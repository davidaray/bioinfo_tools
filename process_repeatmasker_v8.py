#!/usr/bin/env python3
"""
RepeatMasker Overlap Resolution Script
Resolves overlapping repetitive elements using various criteria.
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict
import multiprocessing as mp

class RepeatElement:
    """Represents a single RepeatMasker element"""
    def __init__(self, line: str, line_num: int):
        parts = line.split()
        self.score = int(parts[0])
        self.divergence = float(parts[1])
        self.deletion = float(parts[2])
        self.insertion = float(parts[3])
        self.scaffold = parts[4]
        self.start = int(parts[5])
        self.end = int(parts[6])
        self.left = parts[7]
        self.strand = parts[8]
        self.repeat_name = parts[9]
        self.repeat_class = parts[10]
        self.remaining = ' '.join(parts[11:])
        self.line = line
        self.line_num = line_num
        self.original_line = line
        
        # Handle Simple_repeat class - assign high divergence
        if 'Simple_repeat' in self.repeat_class or 'Low_complexity' in self.repeat_class:
            self.divergence = 999.0
    
    def overlaps(self, other) -> bool:
        """Check if this element overlaps with another"""
        if self.scaffold != other.scaffold:
            return False
        return not (self.end < other.start or self.start > other.end)
    
    def get_overlap_bases(self, other) -> int:
        """Get number of overlapping bases"""
        if not self.overlaps(other):
            return 0
        return min(self.end, other.end) - max(self.start, other.start) + 1
    
    def contains(self, other) -> bool:
        """Check if this element completely contains another"""
        return (self.scaffold == other.scaffold and 
                self.start <= other.start and self.end >= other.end)
    
    def length(self) -> int:
        """Return length of element"""
        return self.end - self.start + 1
    
    def get_class(self) -> str:
        """Extract class from repeat_class field"""
        return self.repeat_class.split('/')[0] if '/' in self.repeat_class else self.repeat_class
    
    def get_family(self) -> str:
        """Extract family from repeat_class field"""
        return self.repeat_class.split('/')[-1] if '/' in self.repeat_class else self.repeat_class
    
    def to_bed(self) -> str:
        """Convert to BED format"""
        strand_bed = '+' if self.strand == '+' else '-'
        # Column 5 is size of hit (end - start + 1)
        hit_size = self.end - self.start + 1
        return f"{self.scaffold}\t{self.start-1}\t{self.end}\t{self.repeat_name}\t{hit_size}\t{strand_bed}\t{self.repeat_class}"
    
    def to_out(self) -> str:
        """Return original RepeatMasker format"""
        return self.original_line


def parse_repeatmasker(filepath: str) -> List[RepeatElement]:
    """Parse RepeatMasker .out file"""
    elements = []
    with open(filepath, 'r') as f:
        for i, line in enumerate(f, 1):
            # Skip first 3 lines (header)
            if i <= 3:
                continue
            line = line.strip()
            if not line:
                continue
            try:
                elements.append(RepeatElement(line, i))
            except (IndexError, ValueError) as e:
                print(f"Warning: Skipping malformed line {i}: {e}", file=sys.stderr)
    return elements


def resolve_overlap_pair(elem1: RepeatElement, elem2: RepeatElement, 
                         resolution: str) -> Tuple[RepeatElement, RepeatElement]:
    """
    Resolve overlap between two elements based on resolution strategy.
    Returns (winner, loser) or (winner, None) if loser is completely removed.
    """
    # Check containment
    e1_contains_e2 = elem1.contains(elem2)
    e2_contains_e1 = elem2.contains(elem1)
    
    if e1_contains_e2 or e2_contains_e1:
        # Containment case
        if resolution == 'higher_score':
            if elem1.score > elem2.score:
                return (elem1, None)
            elif elem2.score > elem1.score:
                return (elem2, None)
            else:
                return (elem1, None)  # First if equal
        
        elif resolution == 'longer_element':
            if e1_contains_e2 and not e2_contains_e1:
                return (elem1, None)
            elif e2_contains_e1 and not e1_contains_e2:
                return (elem2, None)
            else:  # Both equal size
                return (elem1, None)
        
        elif resolution == 'lower_divergence':
            if elem1.divergence < elem2.divergence:
                return (elem1, None)
            elif elem2.divergence < elem1.divergence:
                return (elem2, None)
            else:
                return (elem1, None)  # First if equal
    
    else:
        # Basic overlap - trim the loser
        if resolution == 'higher_score':
            winner = elem1 if elem1.score >= elem2.score else elem2
        elif resolution == 'longer_element':
            winner = elem1 if elem1.length() >= elem2.length() else elem2
        elif resolution == 'lower_divergence':
            winner = elem1 if elem1.divergence <= elem2.divergence else elem2
        
        loser = elem2 if winner == elem1 else elem1
        
        # Trim loser
        if winner.start <= loser.start:
            # Winner is on the left, trim loser's left side
            loser.start = winner.end + 1
        else:
            # Winner is on the right, trim loser's right side
            loser.end = winner.start - 1
        
        # If loser is completely eliminated
        if loser.start > loser.end:
            return (winner, None)
        
        return (winner, loser)


def resolve_overlaps(elements: List[RepeatElement], resolution: str, 
                     progress_callback=None) -> List[RepeatElement]:
    """Resolve all overlaps in list of elements"""
    if not elements:
        return []
    
    # Sort by scaffold and start position
    elements = sorted(elements, key=lambda x: (x.scaffold, x.start, x.line_num))
    
    resolved = []
    i = 0
    total = len(elements)
    
    while i < len(elements):
        if progress_callback and i % 1000 == 0:
            progress_callback(i, total)
        
        current = elements[i]
        j = i + 1
        
        # Find all overlapping elements
        overlapping = []
        while j < len(elements) and elements[j].scaffold == current.scaffold:
            if current.overlaps(elements[j]):
                overlapping.append(j)
                j += 1
            elif elements[j].start > current.end:
                break
            else:
                j += 1
        
        if not overlapping:
            resolved.append(current)
            i += 1
        else:
            # Resolve overlaps iteratively
            group = [current] + [elements[idx] for idx in overlapping]
            resolved_group = resolve_group(group, resolution)
            resolved.extend(resolved_group)
            i = max(overlapping) + 1
    
    return resolved


def resolve_group(group: List[RepeatElement], resolution: str) -> List[RepeatElement]:
    """Resolve overlaps within a group of overlapping elements"""
    if len(group) <= 1:
        return group
    
    # Keep resolving until no overlaps remain
    active = group[:]
    changed = True
    
    while changed:
        changed = False
        new_active = []
        processed = set()
        
        for i, elem1 in enumerate(active):
            if i in processed:
                continue
            
            # Check for overlaps with remaining elements
            found_overlap = False
            for j, elem2 in enumerate(active[i+1:], start=i+1):
                if j in processed:
                    continue
                
                if elem1.overlaps(elem2):
                    winner, loser = resolve_overlap_pair(elem1, elem2, resolution)
                    new_active.append(winner)
                    if loser is not None:
                        # Put loser back for re-evaluation
                        active.append(loser)
                    processed.add(i)
                    processed.add(j)
                    found_overlap = True
                    changed = True
                    break
            
            if not found_overlap:
                new_active.append(elem1)
                processed.add(i)
        
        active = new_active
    
    return active


def filter_elements(elements: List[RepeatElement], dmin: float = None, 
                   dmax: float = None, min_hits: int = None) -> List[RepeatElement]:
    """Filter elements by divergence and minimum hits"""
    # Filter by divergence
    if dmin is not None or dmax is not None:
        filtered = []
        for elem in elements:
            if dmin is not None and elem.divergence < dmin:
                continue
            if dmax is not None and elem.divergence > dmax:
                continue
            filtered.append(elem)
        elements = filtered
    
    # Filter by minimum hits per class/family
    if min_hits is not None:
        class_counts = defaultdict(int)
        for elem in elements:
            class_counts[elem.repeat_class] += 1
        
        elements = [e for e in elements if class_counts[e.repeat_class] >= min_hits]
    
    return elements


def write_output(elements: List[RepeatElement], output_prefix: str, 
                output_type: str, split_by: str = None, progress_dir: str = None):
    """Write output files"""
    
    if progress_dir:
        os.makedirs(progress_dir, exist_ok=True)
        progress_file = os.path.join(progress_dir, f"{output_prefix}_progress.txt")
        with open(progress_file, 'w') as f:
            f.write(f"Writing output files...\n")
    
    if split_by:
        # Group elements by split criterion
        groups = defaultdict(list)
        for elem in elements:
            if split_by == 'class':
                key = elem.get_class()
            elif split_by == 'family':
                key = elem.get_family()
            else:
                key = 'all'
            groups[key].append(elem)
        
        # Write each group
        for key, group_elems in groups.items():
            write_group(group_elems, f"{output_prefix}.{key}", output_type)
    else:
        write_group(elements, output_prefix, output_type)
    
    if progress_dir:
        with open(progress_file, 'a') as f:
            f.write(f"Output complete!\n")


def write_group(elements: List[RepeatElement], prefix: str, output_type: str):
    """Write a group of elements to file(s)"""
    if output_type in ['out', 'both']:
        out_file = f"{prefix}.out"
        with open(out_file, 'w') as f:
            for elem in elements:
                f.write(elem.to_out() + '\n')
        print(f"Wrote {len(elements)} elements to {out_file}")
    
    if output_type in ['bed', 'both']:
        bed_file = f"{prefix}.bed"
        with open(bed_file, 'w') as f:
            for elem in elements:
                f.write(elem.to_bed() + '\n')
        print(f"Wrote {len(elements)} elements to {bed_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Resolve overlaps in RepeatMasker output',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input RepeatMasker .out file')
    parser.add_argument('-s', '--split', choices=['class', 'family'],
                       help='Split output by TE Class or Family')
    parser.add_argument('-ot', '--output-type', required=True,
                       choices=['bed', 'out', 'both'],
                       help='Output file type')
    parser.add_argument('-p', '--prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('-ov', '--overlap-resolution', 
                       choices=['higher_score', 'longer_element', 'lower_divergence'],
                       default='lower_divergence',
                       help='Overlap resolution strategy (default: lower_divergence)')
    parser.add_argument('-m', '--min-hits', type=int,
                       help='Minimum number of hits per TE Class/Family')
    parser.add_argument('-dmax', '--max-divergence', type=float,
                       help='Maximum divergence allowed')
    parser.add_argument('-dmin', '--min-divergence', type=float,
                       help='Minimum divergence allowed')
    parser.add_argument('--progress-dir', 
                       help='Directory for progress updates')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Number of threads (currently single-threaded)')
    
    args = parser.parse_args()
    
    # Set start method for multiprocessing
    mp.set_start_method('spawn', force=True)
    
    print(f"Reading RepeatMasker file: {args.input}")
    elements = parse_repeatmasker(args.input)
    print(f"Loaded {len(elements)} elements")
    
    if args.progress_dir:
        os.makedirs(args.progress_dir, exist_ok=True)
        progress_file = os.path.join(args.progress_dir, f"{args.prefix}_progress.txt")
        with open(progress_file, 'w') as f:
            f.write(f"Loaded {len(elements)} elements\n")
    
    print(f"Filtering elements...")
    elements = filter_elements(elements, args.min_divergence, 
                              args.max_divergence, args.min_hits)
    print(f"After filtering: {len(elements)} elements")
    
    print(f"Resolving overlaps using '{args.overlap_resolution}' strategy...")
    
    def progress(current, total):
        pct = (current / total) * 100
        print(f"Progress: {current}/{total} ({pct:.1f}%)", end='\r')
    
    resolved = resolve_overlaps(elements, args.overlap_resolution, progress)
    print(f"\nResolved: {len(resolved)} elements")
    
    print(f"Writing output files...")
    write_output(resolved, args.prefix, args.output_type, 
                args.split, args.progress_dir)
    
    print("Done!")


if __name__ == '__main__':
    main()

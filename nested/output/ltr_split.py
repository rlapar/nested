#!/usr/bin/env python3

import sys
import glob

def correct_file(fileName):
    output = open("{}_LTR_split.gff".format(fileName[:-4]), 'w')    
    repeat_fragments = []
    ltrs = []
    with open(fileName) as file:
        for _ in range(2):
            output.write(next(file))
        first = file.readline()
        current = first.split("\t")[-1].split(";")[0].split()[-1]
        output.write(first)
        for line in file:
            split_line = line.split("\t")
            chrr = split_line[0]
            idd = split_line[-1].split(";")[0].split()[-1].split("-")[0]
            if idd == current:
                if split_line[2] == "repeat_fragment":
                    repeat_fragments.append([int(split_line[3]), int(split_line[4])])                    
                elif split_line[2] == "long_terminal_repeat":
                    ltrs.append([int(split_line[3]), int(split_line[4])])
                    continue
            else:
                fix_element(output, repeat_fragments, ltrs, chrr, current)
                current = idd                
                repeat_fragments = []
                ltrs = []
            output.write(line)
    fix_element(output, repeat_fragments, ltrs, chrr, current)
    output.close()

def fix_element(output, fragments, ltrs, chrr, current):
    lajna = "{}\tfeature\t{}\t{}\t{}\t.\t.\t.\tID=LTR {} {};Parent=TE_BASE {};name=ltr {}\n"
    ltr = "long_terminal_repeat"
    ltr_n = "ltr_nested_element"
    right, left = ltrs
    left_f = 0
    right_f = 0
    left_nest = 0
    right_nest = 0
    last_frag_end = -1
    for frag in fragments:
        if last_frag_end != -1:
            if frag[0] <= left[1] and within(last_frag_end, left):
                output.write(lajna.format(chrr, ltr_n, last_frag_end + 1, frag[0] - 1, "LEFT", "{}-nest-{}".format(current, left_nest), current, "left"))
                left_nest += 1
            if frag[1] >= right[0] and within(last_frag_end, right):
                output.write(lajna.format(chrr, ltr_n, last_frag_end + 1, frag[0] - 1, "RIGHT", "{}-nest-{}".format(current, right_nest), current, "right"))
                right_nest += 1
        last_frag_end = frag[1]
        
        if intersect(frag, left):
            if left[1] >= frag[1]:
                output.write(lajna.format(chrr, ltr, frag[0], frag[1], "LEFT", "{}-{}".format(current, left_f), current, "left"))
            else:
                output.write(lajna.format(chrr, ltr, frag[0], left[1], "LEFT", "{}-{}".format(current, left_f), current, "left"))
            left_f += 1
        if intersect(frag, right):
            if right[0] >= frag[0]:
                output.write(lajna.format(chrr, ltr, right[0], frag[1], "RIGHT", "{}-{}".format(current, right_f), current, "right"))
            else:
                output.write(lajna.format(chrr, ltr, frag[0], frag[1], "RIGHT", "{}-{}".format(current, right_f), current, "right"))
            right_f += 1

def intersect(a, b):
    return (a[0] >= b[0] and a[0] <= b[1]) or (a[1] >= b[0] and a[1] <= b[1])

def within(point, interval):
    return point >= interval[0] and point <= interval[1]



#for file in sys.argv[1:]:
 #   correct_file(file)

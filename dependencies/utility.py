import inspect
import warnings
import re

def read_fasta(inf, aformat="FIRST", duplicate="replace"):
    data = {}
    with open(inf, "r") as fa:
        name = ""
        for line in fa.readlines():
            if "#" in line:
                continue
            if ">" in line:
                if aformat.upper() == "NCBI":
                    name = re.search(">[a-zA-Z]+_?\d+(\.\d+)*", line).group(0)
                elif aformat.upper() in ["FIRST", "WORD"]:
                    name = line.split()[0]
                else:
                    name = line.strip()
                name = name[1:].strip()
                if name in data.keys():
                    if duplicate.lower() in ["append", "a"]:  # simply add to existing sequence
                        pass
                    elif duplicate.lower() in ["replace", "r"]:  # reset sequence to empty
                        data[name] = ""
                    elif duplicate.lower() in ["separate", "s"]:  # add underscore+number to end of sequence name
                        matches = re.findall("/_\d+$/", name)
                        if matches != None and len(matches) > 0:
                            num = int(max(matches)[1:])
                            name = name[:-len(str(num))] + str(num + 1)
                            data[name] = ""
                        else:
                            name = name + "_2"
                            data[name] = ""
                else:
                    data[name] = ""
            else:
                data[name] = data[name] + line.strip()
    return data

def trim_muts(ntPosnt):
    mut_list = []
    with open(ntPosnt, 'r') as inf:
        for i in inf:
            mut_list.append(i.replace('*',''))

        return [i.replace('\n','') for i in mut_list[1:len(mut_list)]]

def get_mutation_data(ntposnt):
    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    if int(ntposnt[1:-1]) == 1:
        position = int(ntposnt[1:-1])
    else:
        position = int(ntposnt[1:-1]) - 1  # Convert to 0-based index
    return position, (original_nt, mutant_nt)

def convert_position(seq1, seq2, position1, space="-"):
    error = None

    if position1 == 0:
        warnings.warn("\033[93m" + "Position given is not 1-indexed in function " + str(inspect.stack()[1].function) + "\033[0m")
        error = "Position given is not 1-indexed\n" + str(inspect.stack()[1].function)

    i1 = 0
    i2 = 0
    increment = 0
    while i1 < int(position1) and increment < len(seq1) and increment < len(seq2):
        if seq1[increment] != space:
            i1 += 1
        if seq2[increment] != space:
            i2 += 1
        increment += 1
    if not seq1[increment - 1] == seq1.replace(space, "")[i1 - 1]:
        warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq1[increment - 1] + " " + seq1.replace(space, "")[i1 - 1] + "\033[0m")
    elif seq2[increment - 1] != space and not seq2[increment - 1] == seq2.replace(space, "")[i2 - 1]:
        warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq2[increment - 1] + " " + seq2.replace(space, "")[i2 - 1] + "\033[0m")
    if seq2[increment - 1] == space:
        error = "Sequence 1 position aligns with a gap in sequence 2\n" + str(inspect.stack()[1].function)
    return (i2, error)

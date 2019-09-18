#!/usr/bin/python

# Author: Frank Gelhausen

import sys, getopt, subprocess, time
sys.path.append("/usr/local/lib/python2.7/site-packages/")
import RNA
intaRNAPath = "../../IntaRNA_basepairprobs/bin/IntaRNA"

def main(argv):
    seqQ = ""
    seqT = ""
    concQ = 1e-5
    concT = 1e-5

    try:
        opts, args = getopt.getopt(argv, "hq:t:", ["qConc=", "tConc="])
    except getopt.GetoptError:
        print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>]")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>]")
            sys.exit()
        elif opt == "-q":
            seqQ = arg
        elif opt == "-t":
            seqT = arg
        elif opt == "--qConc":
            concQ = float(arg)
        elif opt == "--tConc":
            concT = float(arg)

    (ss, pfT) = RNA.fold_compound(seqT).pf()
    (ss, pfQ) = RNA.fold_compound(seqQ).pf()

    # TODO: why is this required to prevent seg fault?
    (struct, mfe) = RNA.fold(seqT)
    RNA.cvar.cut_point = len(seqT)+1
    (costruct,comfe) = RNA.cofold(seqT + seqQ)
    cmfe = RNA.energy_of_struct(seqT + seqQ, costruct)
    (x,ac,bc,fcab,cf) = RNA.co_pf_fold(seqT + seqQ, struct)
    #print(x, ac, bc, fcab, cf)

    (x,usel1, usel2, fcaa, usel3)= RNA.co_pf_fold(seqT + seqQ, struct)
    RNA.cvar.cut_point = len(seqQ)+1
    (x,usel1, usel2, fcbb, usel3)= RNA.co_pf_fold(seqT + seqQ, struct)

    pfQQ = getPartitionFunction(seqQ, seqQ)
    pfTT = getPartitionFunction(seqT, seqT)
    pfTQ = getPartitionFunction(seqT, seqQ)

    # compare IntaRNA and cofold...
    print(ac, pfT)
    print(bc, pfQ)
    print(fcab, pfTQ)
    print(fcaa, pfTT)
    print(fcbb, pfQQ)

    (AB, AA, BB, A, B) = RNA.get_concentrations(pfTQ, pfTT, pfQQ, pfT, pfQ, concT, concQ)
    AB /= (concT + concQ)
    AA /= (concT + concQ)
    BB /= (concT + concQ)
    A /= (concT + concQ)
    B /= (concT + concQ)

    print(AB, AA, BB, A, B)

def getPartitionFunction(seqT, seqQ):
    res = subprocess.check_output([intaRNAPath, "-q", seqQ, "-t", seqT, "--outMode=E", "--noseed"])
    for line in res.splitlines():
        if "Eall" in line:
            return float(line.split(" ")[1])
    return 0

if __name__ == "__main__":
    main(sys.argv[1:])

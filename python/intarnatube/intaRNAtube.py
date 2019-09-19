#!/usr/bin/python

# Author: Frank Gelhausen

import sys, getopt, subprocess, time
sys.path.append("/usr/local/lib/python2.7/site-packages/")
import RNA
intaRNAPath = "../../../IntaRNA_basepairprobs/bin/IntaRNA"

def main(argv):
    seqT = ""
    seqQ = ""
    cTinit = 1e-5
    cQinit = 1e-5

    try:
        opts, args = getopt.getopt(argv, "hq:t:", ["qConc=", "tConc="])
    except getopt.GetoptError:
        print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>]")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>]")
            sys.exit()
        elif opt == "-t":
            seqT = arg
        elif opt == "-q":
            seqQ = arg
        elif opt == "--tConc":
            cTinit = float(arg)
        elif opt == "--qConc":
            cQinit = float(arg)

    if seqT == "" or seqQ == "":
        print("intaRNAtube.py -q <seq> -t <seq> [--qConc <concentration>] [--tConc <concentration>]")
        sys.exit(2)

    (ss, pfT) = RNA.fold_compound(seqT).pf()
    (ss, pfQ) = RNA.fold_compound(seqQ).pf()

    # TODO: why is this required to prevent seg fault?
    RNA.co_pf_fold("C", None)

    pfTT = getInteractionFeature(seqT, seqT, "Eall")
    pfQQ = getInteractionFeature(seqQ, seqQ, "Eall")
    pfTQ = getInteractionFeature(seqT, seqQ, "Eall")

    (cTQ, cTT, cQQ, cT, cQ) = RNA.get_concentrations(pfTQ, pfTT, pfQQ, pfT, pfQ, cTinit, cQinit)
    cTQ /= (cTinit + cQinit)
    cTT /= (cTinit + cQinit)
    cQQ /= (cTinit + cQinit)
    cT /= (cTinit + cQinit)
    cQ /= (cTinit + cQinit)

    print(cTQ, cTT, cQQ, cT, cQ)

# returns given feature for IntaRNA interaction
def getInteractionFeature(seqT, seqQ, feature):
    res = subprocess.check_output([intaRNAPath, "-q", seqQ, "-t", seqT, "--mode=M", "--outMode=E", "--noseed"])
    for line in res.splitlines():
        if feature in line:
            return float(line.split(" ")[1])
    return 0

if __name__ == "__main__":
    main(sys.argv[1:])

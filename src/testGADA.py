#!/usr/bin/env python3
"""
Examples:
    testGADA.py -i ../data/input.txt -o ../data/output_a0.5T4M5.tsv

Description:
    Test importing eGADA.so to see if it works.
    Each line in the input file is a float number for one probe.
    Output is tab-delimited:
        starting-probe-index, stop-probe-index, no-of-probes, amplitude
"""

import sys
import eGADA


def test_segmentation():
    print("### Testing the C++ eGADA.so module ...\n", flush=True)
    # Pass 1 to eGADA() to enable debugging output.
    # Passing 0 or no passing, i.e. eGADA.eGADA() turns off debugging.
    ins = eGADA.eGADA(1)
    test_vector = [1,1,1,1,0.99,0.99,1,1,0.1,0.1,0.1,0.15]
    segment_ls = ins.run(test_vector, 0.2, 4, 2)
    print(f'\nSegmenting {test_vector} output is:\n \t {segment_ls}.\n', 
            flush=True)
    # check the start, stop, length of the first segment.
    assert segment_ls[0][:3]==[1, 8, 8]
    # check the average amplitude of the first segment.
    assert segment_ls[0][3]==0.9975

    assert segment_ls[1][:3]==[9, 12, 4]
    assert segment_ls[1][3]==0.11249999999999993
    del ins


class call_GADA(object):
    __doc__ = __doc__

    def __init__(self, input_path=None, output_path=None, alpha=0.5, tBackElim=4, \
        minSegLen=5, debug=False):
        # unittest.TestCase.__init__(self,x)
        """
        input_path: Each line in the input file is just one float number for one probe.
        output_path: Path to the output file.
        alpha: alpha in Gamma(alpha; b), the function that controls the prior for the number of breakpoints.
        tBackElim: (amp1-amp2)/stddev in eGADA.
        minSegLen: the minimum number of probes to comprise a segment in eGADA.
        debug: toggle debug mode.
        """
        self.input_path = input_path
        self.output_path = output_path
        self.alpha = alpha
        self.tBackElim = tBackElim
        self.minSegLen = minSegLen
        self.debug = debug


    def run(self):
        """
        """
        #input_path = "/tmp/GADA_ATL8C17938_ATL7C28313_1_5_ATL8C17938_ATL7C28313_1_input"
        inf = open(self.input_path)
        intensity_ls = []
        for line in inf:
            intensity_ls.append(float(line.strip()))
        inf.close()
        sys.stderr.write("%s probes read in.\n"%len(intensity_ls))
        
        ins = eGADA.eGADA()
        segment_ls = ins.run(intensity_ls, self.alpha, self.tBackElim, self.minSegLen)
        del ins
        import pandas as pd
        comment = "# Parameter setting: a=%s,T=%s,minSegLen=%s"%\
            (self.alpha, self.tBackElim, self.minSegLen)
        print(comment)
        header = ['Start', 'Stop', 'Length','Ampl']
        df = pd.DataFrame(data=segment_ls, columns=header)
        df.to_csv(self.output_path, sep='\t', index=False)

if __name__ == '__main__':
    test_segmentation()

    from argparse import ArgumentParser
    ap = ArgumentParser()
                
    ap.add_argument("-i", "--input_path", type=str, required=True,
            help="Each line in the input file is  float number for one probe.")
    ap.add_argument("-o", "--output_path", type=str, required=True,
            help="The path to the output file. Format: "
            "starting-probe-index, stop-probe-index, no-of-probes, amplitude")    
    ap.add_argument("-a", "--alpha", type=float, default=0.5,
            help="Default: %(default)s. "
            "alpha in Gamma(alpha; b), the prior distribution that controls the number of breakpoints.")
    ap.add_argument("-t", "--tBackElim", type=int, default=4,
            help="(amp1-amp2)/stddev in eGADA.")
    ap.add_argument("-m", "--minSegLen", type=int, default=5,
            help='the minimum number of probes in a segment.')
    ap.add_argument("--debug", action='store_true',
            help='Toggle debug mode. Default: %(default)s.')
    args = ap.parse_args()
    instance = call_GADA(args.input_path, args.output_path, \
        alpha=args.alpha, tBackElim=args.tBackElim, \
        minSegLen=args.minSegLen, debug=args.debug)
    instance.run()

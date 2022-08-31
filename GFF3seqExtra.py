## GFF3seqExtra.py ##
## Given geneID List and seq length, identifiy the longest transcripts, output bed file for seq extract ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2022-08-18 ##

## output format 'chrom\tstart\tend\tgenename'

import sys, re

fl_gff3 = sys.argv[1]
fl_geneList = sys.argv[2]
length = int(sys.argv[3])
fl_output = './geneAppend.bed'

dic_longest_transcript = {}
dic_seqRange = {}
## not all gene annotated five_prime_UTR region in gff3 file
dic_seqRange_trans = {}

with open(fl_gff3, 'r') as fh_gff3:
    for line in fh_gff3:
        ## filter comment lines 
        if line[0] != '#':
            line = line.strip('\r\n').split('\t')
            ## construct dic_longest_transcript{}
            if line[2] == 'transcript':
                transInfo = re.split(':|;', line[8])
                gene, transcript = transInfo[3], transInfo[1]
                transLength = int(line[4]) - int(line[3])
                if dic_longest_transcript.get(gene) is None:
                    dic_longest_transcript[gene] = [transcript, transLength]
                elif dic_longest_transcript[gene][1] < transLength:
                    dic_longest_transcript[gene] = [transcript, transLength]
                else:
                    pass
                ## constract dic_seqRange_trans{}
                direction = line[6]
                if direction == '-':
                    chrom, start, end = line[0], int(line[4]), int(line[4]) + length
                else:
                    chrom, start, end = line[0], int(line[3]), (int(line[3]) - length, 0)[int(line[3]) - length < 0]
                dic_seqRange_trans[transcript] = [chrom, str(min(start, end)), str(max(start, end))]

            ## construct dic_seqRange{}
            if line[2] == 'five_prime_UTR':
                transcript = re.split(':', line[8])[1]
                direction = line[6]
                if direction == '-':
                    chrom, start, end = line[0], int(line[3]), int(line[3]) + length
                else:
                    chrom, start, end = line[0], int(line[4]), (int(line[4]) - length, 0)[int(line[4]) - length < 0]
                dic_seqRange[transcript] = [chrom, str(min(start, end)), str(max(start, end))]

with open(fl_geneList, 'r') as fh_geneList:
    res = []
    for line in fh_geneList:
        gene = line.strip('\r\n')
        if dic_longest_transcript.get(gene) is not None:
            transcript = dic_longest_transcript[gene][0]
            if dic_seqRange.get(transcript) is not None:
                transcriptInfo = dic_seqRange[transcript]
                outputfmt = transcriptInfo + [gene]
                res.append(outputfmt)
            else:
                transcriptInfo = dic_seqRange_trans[transcript]
                outputfmt = transcriptInfo + [gene]
                res.append(outputfmt)

with open(fl_output, 'w') as fh_output:
    res = ['\t'.join(i) for i in res]
    fh_output.write('\n'.join(res))

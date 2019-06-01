#!/usr/bin/env python

#Using Dovetail table, convert coordinates in bedfile to new genome coordinates
from optparse import OptionParser
import gzip

def parse_options():
    USAGE = """
    convert_coord_gff.py --coordinate table
                --gff input
                --out
    """
    parser = OptionParser(USAGE)

    parser.add_option('--dovetail_table', default=None, help='path to table of coordinates from Dovetail Assembly')
    parser.add_option('--gff_file', default=None, help='path to the gff file of annotations -expects 9 column format')
    parser.add_option('--output_gff', default=None, help='path to output file')

    (options,args)=parser.parse_args()
    return(options)

    
def change_coords(origA, out):
    
    for rec in origA:
        if rec.startswith('##sequence'):
            continue
        elif rec.startswith('#'):
            out.write(rec)
            continue
        else:
            recLine = rec.strip().split("\t")
            BdSc = recLine[0]
            BdSt = int(recLine[3])
            BdEd = int(recLine[4])
            BdStrd = recLine[6]
        
            if options.dovetail_table.endswith('gz'):
                dovetail = gzip.open(options.dovetail_table)
            else:
                dovetail = open(options.dovetail_table)    
            dovetail.readline()
            dovetail.readline() #skip the header
            for line in dovetail:
                myLine = line.strip().split()
                DTsc = myLine[0]
                ORsc = myLine[1]
                ORst = int(myLine[2])
                OREd = int(myLine[3])
                strand = myLine[4]
                DT_start = int(myLine[5])
                DT_end = int(myLine[6])
                
                NewSC = "NA"
                NewSt = "NA"
                NewEd = "NA"
                New_strand = "NA"
                AdSt = "NA"
                AdEd = "NA"
             
                if BdSc != ORsc: continue
                
                else:    
                    #print NewSc, "\t", BdSc,"\n"
                    if BdSt >= ORst and BdEd <= OREd:
                        NewSc = DTsc
                        #print NewSc, "\t", BdSc, "\t", BdSt, "\t", ORst, "\n"
                        if ORst > 0:
                            AdSt = int(BdSt) - int(ORst)
                            AdEd = int(OREd) - int(BdEd)
                            if strand == "-" and BdStrd == "+":
                                NewSt = (int(DT_start) + int(AdEd)) + 1
                                #print NewSt, "\n
                                NewEd = (int(DT_end) - int(AdSt)) + 1
                                New_strand = "-"
                                out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                                out.write('\t'.join(recLine[7:]) + '\n')
                                break
                            elif strand == "-" and BdStrd == "-":
                                NewSt = (int(DT_start) + int(AdEd)) + 1
                                #print NewSt, "\n"
                                NewEd = (int(DT_end) - int(AdSt)) + 1
                                New_strand = "+"
                                out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                                out.write('\t'.join(recLine[7:]) + '\n')
                                break
                            elif strand == "+" and BdStrd == "-":
                                NewSt = int(DT_start) + int(AdSt)
                                NewEd = int(DT_end) - int(AdEd)
                                New_strand = "-"
                                out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                                out.write('\t'.join(recLine[7:]) + '\n')
                                break
                            elif strand == "+" and BdStrd == "+":
                                NewSt = int(DT_start) + int(AdSt)
                                NewEd = int(DT_end) - int(AdEd)
                                New_strand = "+"
                                out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                                out.write('\t'.join(recLine[7:]) + '\n')
                                break
                            elif strand == "-" and BdStrd == ".":
                                NewSt = (int(DT_start) + int(AdEd)) + 1
                                NewEd = (int(DT_end) - int(AdSt)) + 1
                                New_strand = "-"
                                out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                                out.write('\t'.join(recLine[7:]) + '\n')
                                break
                            elif strand == "+" and BdStrd == ".":
                                NewSt = int(DT_start) + int(AdSt)
                                NewEd = int(DT_end) - int(AdEd)
                                New_strand = "+"
                                out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                                out.write('\t'.join(recLine[7:]) + '\n')
                                break
                            else:
                                print "No matches for", BdSc, recLine[1], "\t", recLine[2], "\t", BdSt, "\t", BdEd, "\t", recLine[5], "\t", BdStrd, "\t", recLine[7:], "!"
                                break
                            #print NewSc, "\t", BdSc, "\t", BdSt, "\t", ORst,"\t", BdEd, "\t", OREd, "\n"
                            #print '\t'.join(recLine), "\n"
                        elif strand == "-" and BdStrd == "+":
                            NewSt = (int(DT_end) - int(BdEd)) + 1
                            #print NewSt, "\n"
                            NewEd = (int(DT_end) - int(BdSt)) + 1
                            #print NewEd, "\n"
                            New_strand = "-"
                            out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                            out.write('\t'.join(recLine[7:]) + '\n')
                            break
                        elif strand == "-" and BdStrd == "-":
                            NewSt= (int(DT_end) - int(BdEd)) + 1
                            #print NewSt, "\n"
                            NewEd = (int(DT_end) - int(BdSt)) + 1
                            #print NewEd, "\n"
                            New_strand = "+"
                            out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                            out.write('\t'.join(recLine[7:]) + '\n')
                            break
                        elif strand == "+" and BdStrd == "-":
                            NewSt = int(DT_start) + int(BdSt)
                            NewEd = int(DT_start) + int(BdEd)
                            New_strand = "-"
                            out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                            out.write('\t'.join(recLine[7:]) + '\n')
                            break
                        elif strand == "+" and BdStrd == "+":
                            NewSt = int(DT_start) + int(BdSt)
                            NewEd = int(DT_start) + int(BdEd)
                            New_strand = "+"
                            out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                            out.write('\t'.join(recLine[7:]) + '\n')
                            break
                        elif strand == "-" and BdStrd == ".":
                            NewSt = (int(DT_end) - int(BdEd)) + 1
                            #print NewSt, "\n"
                            NewEd = (int(DT_end) - int(BdSt)) + 1
                            #print NewEd, "\n"
                            New_strand = "-"
                            out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                            out.write('\t'.join(recLine[7:]) + '\n')
                            break
                        elif strand == "+" and BdStrd == ".":
                            NewSt = int(DT_start) + int(BdSt)
                            NewEd = int(DT_start) + int(BdEd)
                            New_strand = "+"
                            out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                            out.write('\t'.join(recLine[7:]) + '\n')
                            break
                        else:
                            print "No matches for", BdSc, recLine[1], "\t", recLine[2], "\t", BdSt, "\t", BdEd, "\t", recLine[5], "\t", BdStrd, "\t", recLine[7:], "!"
                            break
                    else: continue
            
                        #print NewSt
                        #print NewEd
                        #print New_strand
            
            
            dovetail.close()
        
'''if BdSc == ORsc and BdSt >= ORst and BdEd <= OREd:
                    NewSc = DTsc
                    if strand == "-" and BdStrd == "+":
                        NewSt = (DT_end - BdEd)
                        NewEd = (DT_end - BdSt)
                        New_strand = "-"
                    elif strand == "-" and BdStrd == "-":
                        NewSt = (DT_end - BdEd)
                        NewEd = (DT_end - BdSt)
                        New_strand = "+"
                    elif strand == "+" and BdStrd == "-":
                        NewSt = (DT_start + BdSt)
                        NewEd = (DT_start + BdEd)
                        New_strand = "-"
                    elif strand == "+" and BdStrd == "+":
                        NewSt = (DT_start + BdSt)
                        NewEd = (DT_start + BdEd)
                        New_strand = "+"
                    out.write('%s\t%s\t%s\t%i\t%i\t%s\t%s\t' %(NewSc, recLine[1], recLine[2], NewSt, NewEd, recLine[5], New_strand))
                    out.write('\t'.join(recLine[7:]) + '\n')
            else:
                continue '''
                    
if __name__ == '__main__':
  #load parameters and files
    options = parse_options()
    if options.gff_file.endswith('gz'):
        gff_file = gzip.open(options.gff_file)
    else:
        gff_file = open(options.gff_file)
        
    outfile = open(options.output_gff + '.gff', 'w')
    
    change_coords(gff_file, outfile)
    
    gff_file.close()
    outfile.close()
    
                    
            
            
            
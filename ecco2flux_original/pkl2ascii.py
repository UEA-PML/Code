#!/usr/bin/env python
import sys
import pandas as pd
import argparse

# Description:
#    Python code to convert a level0 pkl format dataframe into an output ASCII file
# USAGE: pkl2ascii.py -L0 <L0 input filename> -o <output filename>
#
def main():

   parser = argparse.ArgumentParser(
      description=__doc__,
      formatter_class=argparse.RawDescriptionHelpFormatter
   )

   parser.add_argument("-o", "--ofile", type=str, default='', help="Output filename")
   parser.add_argument("-L0", "--L0", type=str, default='', help="L0 input filestem")
   args = parser.parse_args()

   if not(args.L0 or args.ofile):
      print "USAGE: pkl2ascii.py -L0 <L0 input filename> -o <output filename>"
      sys.exit()
   
   if args.L0:
      ifilename = args.L0
      print "Input filename: ", ifilename
      freq_data_matrix = pd.read_pickle(ifilename)

   if args.ofile:   
      ofilename = args.ofile
      print " Output filename: ", ofilename
      FILE = open(ofilename, 'w')
      header_str = ['u(m/s)','v(m/s)','w(m/s)','T(C)',
                    'rotx(deg)','roty(deg)','rotz(deg)','ax(m/s2)','ay(m/s2)','az(m/s2)',
                    'Pic_CH4_conc(ppm)','Pic_CO2_conc(ppm)','Switch','Licor_CO2_dens','Licor_H20_dens',
                    'Licor_T(degC)','Licor_P(mb)','Licor_diag','Latititude(deg)','Longitude(deg)',
                    'Heading(deg)','Speed(knots)']

      FILE.write("%s\n" % (freq_data_matrix.to_string(header=header_str)));
      FILE.close()


   sys.exit()
      
if __name__=='__main__':
   main()

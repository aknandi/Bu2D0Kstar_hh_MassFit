# Calculate R(D*K*/DK*) and frac010 for the binned fit
import sys

# inputs
data={}
inputrange=sys.argv[1]
tag=sys.argv[2]
dir = '../Settings/PDFShapes/Gen/'+tag+'/'+inputrange

# output files
output_file_ratios = open(dir+'/ratiosFIXED_'+inputrange+'_corr.txt','w')
output_file_partreco = open(dir+'/partrecoFIXED_corr.txt','w')

# input files
loc_partrecoFIXED = dir+'/partrecoFIXED.txt'
loc_gentotals = '../output/GenTotals_binnedfit_'+inputrange+'.txt'

# Get the info we need and store in 'data'
input_file = open(loc_partrecoFIXED,'r')
for line in input_file:
  entry = line.split()
  if entry[0]=='frac010_bs_'+inputrange:
    data[entry[0]]=float(entry[1])
input_file.close()

#input_file = open('../Settings/PDFShapes/Gen/'+inputrange+'/ratiosFIXED_'+inputrange+'.txt','r')
#for line in input_file:
#  entry = line.split()
#  if 'ratio_bs_dstkst_both_mix_merge'==entry[0]:
#    data[entry[0]]=float(entry[1])
#  if 'ratio_bs_drho_both_mix_merge'==entry[0]:
#    data[entry[0]]=float(entry[1])
#input_file.close()
#

#for item in data:
#  print item, data[item]

# the following numbers obtained by running the fit and printing out the calculated values of G, P = BR * Eff * IntegFrac
data['G001_4900']=0.00541477
data['P001_4900']=0.0089167
data['G010_4900']=0.00875386
data['P010_4900']=0.0142729
data['G001_5160']=0.00311149
data['P001_5160']=0.00423723
data['G010_5160']=0.0047879
data['P010_5160']=0.00732387
data['G001_5200']=0.00256833
data['P001_5200']=0.000652012
data['G010_5200']=0.00314023
data['P010_5200']=0.00278808

# Calculate the new alpha010
fAlpha = data['frac010_bs_'+inputrange] / (1-data['frac010_bs_'+inputrange])
f4900 = (data['G001_'+inputrange] + data['P001_'+inputrange]) / (data['G010_'+inputrange] + data['P010_'+inputrange]) 
f5200 = (data['G010_5200'] + data['P010_5200']) / (data['G001_5200'] + data['P001_5200']) 
frac = fAlpha * f4900 * f5200
alpha2 = frac / (1+frac)

# Calculate the new R(D*K*/DK*)
#alpha1 = data['frac010_bs_'+inputrange] 
#
#up = alpha2 * (data['G010_5200'] + data['P010_5200']) + (1-alpha2) * (data['G001_5200'] + data['P001_5200'])
#down = alpha1 * (data['G010_'+inputrange] + data['P010_'+inputrange]) + (1-alpha1) * (data['G001_'+inputrange] + data['P001_'+inputrange])
#ratio = data['ratio_bs_dstkst_both_mix_merge']
#ratio2 = ratio * up / down
#

# Calculate the new R(D*K*/DK*) just based on the yield outputs
input_file = open(loc_gentotals,'r')
for line in input_file:
  entry = line.split()
  if 'N_bs_d2kspipi_both_mix'==entry[0]:
    n_bs = float(entry[1])
  if 'N_bs_dstkst_d2kspipi_both_mix'==entry[0]:
    n_dstkst = float(entry[1])
  if 'N_drho_d2kspipi_both_mix'==entry[0]:
    n_drho=float(entry[1])
input_file.close()

ratio_dstkst = n_dstkst / n_bs
ratio_drho = n_drho / n_bs

# Output to necessary text files
#print 'ratio_bs_dstkst_both_mix_merge {:.4}'.format(ratio2)

output_file_ratios.write('ratio_bs_dstkst_both_mix_merge {:.4}\n'.format(ratio_dstkst))
output_file_ratios.write('ratio_bs_dstkst_both_mix_merge_LimL {:.4}\n'.format(ratio_dstkst))
output_file_ratios.write('ratio_bs_dstkst_both_mix_merge_LimU {:.4}\n'.format(ratio_dstkst))
output_file_ratios.write('ratio_bs_drho_both_mix_merge {:.4}\n'.format(ratio_drho))
output_file_ratios.write('ratio_bs_drho_both_mix_merge_LimL {:.4}\n'.format(ratio_drho))
output_file_ratios.write('ratio_bs_drho_both_mix_merge_LimU {:.4}\n'.format(ratio_drho))
output_file_ratios.close()

#print 'frac010_bs_5200 {:.4}'.format(alpha2)
input_file = open(loc_partrecoFIXED,'r')
for line in input_file:
  if not 'frac010_bs' in line:
    output_file_partreco.write(line)
# Note the following is always 5200 because it is for the binned fit
output_file_partreco.write('frac010_bs_5200 {:.4}\n'.format(alpha2))
output_file_partreco.write('frac010_bs_5200_LimL {:.4}\n'.format(alpha2))
output_file_partreco.write('frac010_bs_5200_LimU {:.4}\n'.format(alpha2))
output_file_partreco.close()


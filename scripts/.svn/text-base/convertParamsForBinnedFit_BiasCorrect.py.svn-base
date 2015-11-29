# Calculate R(D*K*/DK*) and frac010 for the binned fit
import sys

# inputs
data={}
vars={}
inputrange=sys.argv[1]

output_file = open('../output/values_biascorr_rangecorr_'+inputrange+'.txt','w')

# Get the info we need and store in 'data'
input_file = open('../output/values_biascorr_'+inputrange+'.txt','r')
for line in input_file:
  entry = line.split()
  vars[entry[0]]=float(entry[1])
input_file.close()

# Get the range corrected values of ratios
ratio_dstkst_rangecorr=0
ratio_drho_rangecorr=0
input_file = open('../Settings/PDFShapes/Gen/'+inputrange+'/ratiosFIXED_'+inputrange+'_corr.txt','r')
for line in input_file:
  entry = line.split()
  if(entry[0]=='ratio_bs_dstkst_both_mix_merge'):
    ratio_dstkst_rangecorr=float(entry[1])
  if(entry[0]=='ratio_bs_drho_both_mix_merge'):
    ratio_drho_rangecorr=float(entry[1])
input_file.close()

# Get the values of ratios
ratio_dstkst=0
ratio_drho=0
input_file = open('../Settings/PDFShapes/Gen/'+inputrange+'/ratiosFIXED_'+inputrange+'.txt','r')
for line in input_file:
  entry = line.split()
  if(entry[0]=='ratio_bs_dstkst_both_mix_merge'):
    ratio_dstkst=float(entry[1])
  if(entry[0]=='ratio_bs_drho_both_mix_merge'):
    ratio_drho=float(entry[1])
input_file.close()

# Get the yields
n_bs_rangecorr, n_bd_rangecorr, n_comb_rangecorr = 0, 0, 0
input_file = open('../output/GenTotals_binnedfit_'+inputrange+'.txt','r')
for line in input_file:
  entry = line.split()
  if(entry[0]=='N_bs_d2kspipi_both_mix'):
    n_bs_rangecorr=float(entry[1])
  if(entry[0]=='N_bd_d2kspipi_both_mix'):
    n_bd_rangecorr=float(entry[1])
  if(entry[0]=='N_comb_d2kspipi_both_mix'):
    n_comb_rangecorr=float(entry[1])
input_file.close()

# Get the reduced yields to 5200
n_bs, n_bd, n_comb =0, 0, 0
input_file = open('../output/GenTotals_'+inputrange+'.txt','r')
for line in input_file:
  entry = line.split()
  if(entry[0]=='N_bs_d2kspipi_both_mix'):
    n_bs=float(entry[1])
  if(entry[0]=='N_bd_d2kspipi_both_mix'):
    n_bd=float(entry[1])
  if(entry[0]=='N_comb_d2kspipi_both_mix'):
    n_comb=float(entry[1])
input_file.close()


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
fAlpha = vars['frac010_bs_'+inputrange] / (1-vars['frac010_bs_'+inputrange])
f4900 = (data['G001_'+inputrange] + data['P001_'+inputrange]) / (data['G010_'+inputrange] + data['P010_'+inputrange]) 
f5200 = (data['G010_5200'] + data['P010_5200']) / (data['G001_5200'] + data['P001_5200']) 
frac = fAlpha * f4900 * f5200
alpha2 = frac / (1+frac)

# Calculate the new ratios
newratio_dstkst = vars['ratio_bs_dstkst_both_mix_merge'] / ratio_dstkst * ratio_dstkst_rangecorr
newratio_drho = vars['ratio_bs_drho_both_mix_merge'] / ratio_drho * ratio_drho_rangecorr
new_n_bd = vars['n_bd_fit_d2kspipi_both_mix_merge'] / n_bd * n_bd_rangecorr
new_n_bs = vars['n_bs_fit_d2kspipi_both_mix_merge'] / n_bs * n_bs_rangecorr
new_n_comb = vars['n_comb_d2kspipi_both_mix_merge'] / n_comb * n_comb_rangecorr

# put back into the array
vars['frac010_bs_'+inputrange] = alpha2
vars['ratio_bs_dstkst_both_mix_merge'] = newratio_dstkst
vars['ratio_bs_drho_both_mix_merge'] = newratio_drho
vars['n_bs_fit_d2kspipi_both_mix_merge'] = new_n_bs
vars['n_bd_fit_d2kspipi_both_mix_merge'] = new_n_bd
vars['n_comb_d2kspipi_both_mix_merge'] = new_n_comb

# Output to necessary text file
for entry in vars:
  output_file.write(entry+' {:.6f}\n'.format(vars[entry]))
output_file.close()


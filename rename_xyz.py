import os

pstring = os.path.join('xyz_files', 'DMPr_CG_TABT_XYZ')

for f in os.listdir(pstring):
    os.rename(os.path.join(pstring, f), 
    os.path.join(pstring, f.replace('1,3-DMP', 'DMPr')))
                 
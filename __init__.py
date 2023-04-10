# 2021 01 20
#import feedparser
from flask import Flask,g,request, Response,make_response, render_template, jsonify, redirect, url_for
from torch import FloatStorage
from flask_dropzone import Dropzone
#from werkzeug import secure_filename
from werkzeug.utils import secure_filename
import os
from subprocess import *
from mpl_toolkits.mplot3d import Axes3D

import matplotlib
matplotlib.use('Agg') # no UI backend
import matplotlib.pyplot as plt

import numpy as np
import math
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

from pymongo import MongoClient
from pymongo.cursor import CursorType

from glob import glob

from werkzeug.datastructures import ImmutableMultiDict
import time
import make_html
import Field_names

import bcrypt 
from vasprun import vasprun

#------------------------------------------------------------------------------------------------
def remove_side(ori, Left, Right ):
    item_string = str(ori)
    item_string = item_string.lstrip(Left)
    item_string = item_string.rstrip(Right)
    item_string = item_string.strip()

    return item_string

#------------------------------------------------------------------------------------------------
def replace_a_line(f_name, n, new_string):
    file0 = open(f_name, 'r')
    f_temp= f_name + "_temp"
    file1 = open(f_temp, 'w')
    i=0
    while True:         
         line = file0.readline() 
         if i==n: line = new_string + "\n"
         if not line : break
         else : file1.write(line)  
         i += 1
    file0.close()
    file1.close()
    command_string  ="rm {0}; mv {1} {2}".format(f_name, f_temp, f_name); os.system(command_string)
#------------------------------------------------------------------------------------------------
def graph_process(t_dir):    
    
    graph_f_name  = t_dir + "/tDOS.png"

    #display_graph = "display " + graph_f_name 
    plt.savefig(graph_f_name)
    #os.system(display_graph)

# Read_xml -----------------------------------------------------------------
def Read_xml(t_dir, c_mode):
    global DOS_title, N_p_title, Band_title
    global Nk_per_div, N_bands
    global Band_max, Band_min 
    global k_acc
    global N_Div
    global E_Fermi, V_crystal
    global basis
    global Is_pDOS
    global N_atoms
    global N_E

    global Is_Band
    global Is_spin_dn
    


    DOS_title = []
    Band_title = []

    DOS_ON = False
    
    tDOS_ON = False
    tDOS_up_ON = False
    tDOS_dn_ON = False
    
    pDOS_ON = False    
    pDOS_up_ON = False
    pDOS_dn_ON = False

    kPOI_ON = False
    weig_ON = False
    Band_ON = False 
    parr_up_ON = False
    parr_dn_ON = False

    Barr_ON = False
    Rdiv_ON = False
    Proj_ON = False
    pBand_ON = False
    pBarr_ON = False
    FinalPos_ON = False
    basis_ON = False
    N_atoms_ON = False
    kPOI_read = False

    Is_pDOS = False
    Is_Band = False
    Is_pBand = False
    Is_spin_dn = False
   
    EeV  = []
    tDOS_up = []; tDOS_dn = []; pDOS_up = []; pDOS_dn = []
    k    = []; kx   = []; ky   = []; kz   = []; wei  = []; k_acc = []
    Band = []
    pBand = []
    pBand_title = []
    basis = []
    basis = [[] for _ in range(3)]
    
    filename_string  = t_dir + "/vasprun.xml"; file0 = open(filename_string, 'r')   
    
    N_E = 0
    N_p_title = 0
    N_pBand_title = 0
    N_kP = 0

    Nk_per_div =-1

    N_bands = 0

    Band_i =0
    N_Div =0

    while True:
       
        line = file0.readline().strip()

        if not line : 
            break        
        # Read pBands
        elif FinalPos_ON == True and '"volume"' in line: line = line.replace('<i name="volume">','').replace('</i>','').strip(); V_crystal = float(line) ; FinalPos_ON   = False; basis_ON = False
        elif basis_ON    == True and      '<v>' in line: 
            line = line.replace('<v>','').replace('</v>','').strip(); 
            line = ' '.join(line.split()); line_list = line.split(' '); 
            basis[basis_i].append(np.double(line_list[0])); 
            basis[basis_i].append(np.double(line_list[1])); 
            basis[basis_i].append(np.double(line_list[2])); 
            basis_i +=1  
        elif FinalPos_ON == True and  '"basis"' in line: basis_ON   = True

        elif '"finalpos"' in line : FinalPos_ON   = True ; basis_i=0

        elif Proj_ON == True and line=="</projected>"   : Proj_ON   = False; pBand_ON   = False

        elif pBand_ON == True and '<r>' in line         : 
            if pBarr_ON == False: pBand=  [[[] for _ in range(N_bands)]  for _ in range(N_pBand_title)]; pBarr_ON = True 
            line = ' '.join(line.split()); line_list = line.split(' ')
            for i in range(N_pBand_title): pBand[i][Band_i].append(float(line_list[i+1]))
            Band_i += 1
            if Band_i == N_bands: Band_i = 0  
        elif Proj_ON == True and Band_ON == False and line=="<array>" : pBand_ON   = True ; Band_i = 0       
        elif Proj_ON == True and Band_ON == False and '<field>' in line : line = line.replace('<field>','').replace('</field>','').strip(); pBand_title.append(line); N_pBand_title +=1      
        
        # Read Bands
        elif Band_ON == True and line=="</eigenvalues>" : Band_ON   = False
        elif Band_ON == True and '<r>' in line          : 
            if Barr_ON == False: Band=  [[] for _ in range(N_bands)]; Barr_ON = True 
            line = ' '.join(line.split()); line_list = line.split(' '); Band[Band_i].append(float(line_list[1])); Band_i += 1
            if Band_i == N_bands: Band_i = 0  
        elif line    =='<eigenvalues>'                  : Band_ON   = True; Band_i = 0  ; Is_Band = True

        # Read on Band and pBand
        elif line    =='<projected>'                    : Proj_ON   = True      ; Is_pBand = True
          
        # Read partial DOS
        elif DOS_ON == True and line =="</dos>"              : DOS_ON = False

        elif pDOS_ON   == True and line    =="</partial>"   : pDOS_ON   = False          
        elif pDOS_dn_ON == True and line=="</set>"         : pDOS_dn_ON   = False
        elif pDOS_dn_ON == True and     '<r>' in line      : 
            line = ' '.join(line.split()); line_list = line.split(' ')
            if parr_dn_ON == False: pDOS_dn=  [[] for _ in range(N_p_title-1)]; parr_dn_ON = True 
            for d_i in range(N_p_title-1): pDOS_dn[d_i].append(line_list[d_i+2]) 
        elif pDOS_ON == True and '"spin 2"' in line     : pDOS_dn_ON = True


        elif pDOS_up_ON == True and line=="</set>"         : pDOS_up_ON   = False
        elif pDOS_up_ON == True and     '<r>' in line      : 
            line = ' '.join(line.split()); line_list = line.split(' ')
            if parr_up_ON == False: pDOS_up=  [[] for _ in range(N_p_title-1)]; parr_up_ON = True 
            for d_i in range(N_p_title-1): pDOS_up[d_i].append(line_list[d_i+2]) 
        elif pDOS_ON == True and '"spin 1"' in line     : pDOS_up_ON = True

        elif pDOS_ON == True and '<field>' in line      : line = line.replace('<field>','').replace('</field>','').strip(); DOS_title.append(line); N_p_title +=1
        elif line    =="<partial>"                      : pDOS_ON   = True ; Is_pDOS = True
        # Read total DOS
        
        elif tDOS_ON    == True and line=="</total>"       : tDOS_ON      = False
        elif tDOS_dn_ON == True and line=="</set>"         : tDOS_dn_ON   = False

        elif tDOS_dn_ON == True and '<r>' in line          : line = ' '.join(line.split()); line_list = line.split(' '); tDOS_dn.append(line_list[2]); 
        elif tDOS_ON    == True and '"spin 2"' in line     : tDOS_dn_ON = True; Is_spin_dn = True

        elif tDOS_up_ON == True and line=="</set>"         : tDOS_up_ON = False
        elif tDOS_up_ON == True and '<r>' in line          : line = ' '.join(line.split()); line_list = line.split(' '); EeV.append(float(line_list[1])); tDOS_up.append(line_list[2]); N_E+=1

        elif tDOS_ON == True and '"spin 1"' in line        : tDOS_up_ON = True
        elif DOS_ON  == True and line == '<total>'         : tDOS_ON = True 
        elif DOS_ON  == True and 'efermi' in line          : line = line.replace('<i name="efermi">','').replace('</i>','').strip(); E_Fermi = float(line) 
        elif line =="<dos>"                                : DOS_ON = True 
      
        # Read N_bands
        elif '"NBANDS"' in line                         : line = line.replace('<i type="int" name="NBANDS">','').replace('</i>','').strip(); N_bands = int(line)  
        
        # Read weights
        elif weig_ON == True and line=="</varray>"      : weig_ON   = False
        elif weig_ON == True and '<v>' in line          : line = ' '.join(line.split()); line_list = line.split(' '); wei.append(line_list[1])  
        elif line    =='<varray name="weights" >'       : weig_ON   = True
        
        # Read k points
        elif kPOI_ON == True and line=="</varray>"      : kPOI_ON   = False; kPOI_read = True
        elif kPOI_ON == True and '<v>' in line          : 
            line = ' '.join(line.split()); 
            line_list = line.split(' '); 
            kx.append(float(line_list[1])); 
            ky.append(float(line_list[2])); 
            kz.append(float(line_list[3])); 
            #k.append(math.sqrt(kx[N_kP]**2+ky[N_kP]**2+kz[N_kP]**2)); 
            N_kP +=1 
        elif line=='<varray name="kpointlist" >' and kPOI_read == False : kPOI_ON   = True

        elif Rdiv_ON == True and line=="</generation>"  : Rdiv_ON = False
        elif Rdiv_ON == True and '<v>' in line          : N_Div += 1        
        elif '<i name="divisions" type="int">' in line  : print('line'+line); line = line.replace('<i name="divisions" type="int">','').replace('</i>','').strip(); line_list = line.split(' '); Nk_per_div = int(line_list[0]); Rdiv_ON= True

        elif N_atoms_ON==True and '<atoms>' in line     : line = line.replace('<atoms>','').replace('</atoms>','').strip()  ; N_atoms = int(line); N_atoms_ON = False
        elif line=='<atominfo>'        : N_atoms_ON   = True
   
    file0.close()
  
    #  save DOS ----------------------------------------------------------------------------- 
    if c_mode ==1:
        
        filename_string  = t_dir + "/tDOS.dat"   ; file0 = open(filename_string, 'w')
        write_string = "energy,total_up"
        if Is_spin_dn: write_string += ",total_dn"

        write_string +='\n'; 

        for i in range(N_E)           : 
            EeV[i] -= E_Fermi
            write_string  += str(EeV[i]) + ',' + str(tDOS_up[i]) 
            if Is_spin_dn : write_string  += ',' + str(tDOS_dn[i]) 
            write_string  += '\n'            
            
        file0.write(write_string)
        file0.close()

        if Is_pDOS :
            filename_string  = t_dir + "/pDOS.dat"   ; file0 = open(filename_string, 'w')
            
            write_string = DOS_title[0]
            for i in range(1, N_p_title) :
                write_string += ","+ DOS_title[i] + "_up"
                if Is_spin_dn : write_string += ","+ DOS_title[i] + "_dn"
            write_string +='\n'

            for i in range(N_E) :             
                write_string  += str(EeV[i])
                for j in range(N_p_title-1) : 
                    write_string += ',' + str(pDOS_up[j][i])        
                    if Is_spin_dn : write_string += ',' + str(pDOS_dn[j][i])
                write_string += '\n'
        
            file0.write(write_string)
                
        file0.close()



    if c_mode ==3:        
        filename_string  = t_dir + "/tDOS.dat"   ; file0 = open(filename_string, 'w')
        write_string = "energy,total\n"
        file0.write(write_string)

        for i in range(N_E)           : 
            EeV[i] -= E_Fermi
            write_string  = str(EeV[i])+','+ str(tDOS_up[i])            
            write_string += '\n'
            file0.write(write_string)
        file0.close()

    
    #  save Bands ----------------------------------------------------------------------------- 
    if Is_Band == True:    
        filename_string  = t_dir + "/Bands.dat" ; file0 = open(filename_string, 'w')
        write_string = 'k'    
        for i in range(N_bands): Band_title.append('Band'+ str(i)); write_string += ','+ Band_title[i] 
        write_string += '\n' 
        file0.write(write_string)

        k_acc_v =0

        Band_max = -1.0e100
        Band_min =  1.0e100

        print('N_kP',N_kP)
        print('N_bands',N_bands)

        for i in range(N_kP) :
            if i==0 or i%Nk_per_div ==0 : pass
            else                        : k_acc_v += math.sqrt((kx[i]-kx[i-1])**2+(ky[i]-ky[i-1])**2+(kz[i]-kz[i-1])**2);                     

            k_acc.append(k_acc_v)           
            write_string  = str(k_acc_v) 
            for j in range(N_bands): 
                #print(i,j,Band[j][i])
                Band[j][i] -= E_Fermi
                write_string += ','+ str(Band[j][i]); 
                if Band[j][i]> Band_max: Band_max = Band[j][i]
                if Band[j][i]< Band_min: Band_min = Band[j][i] 
            write_string +='\n' 
            file0.write(write_string)
        file0.close()

    #  save pBands ----------------------------------------------------------------------------- 
    if Is_pBand == True:          
        for k in range(N_pBand_title):    
            filename_string  = t_dir + "/pBands_" + str(k) + ".dat" ; file0 = open(filename_string, 'w')   
            write_string = 'k'    
            for i in range(N_bands): Band_title.append('Band'+ str(i)); write_string += ','+ Band_title[i]    
            write_string += '\n'     
            file0.write(write_string)
            #Band_max = -1.0e100
            #Band_min =  1.0e100
            
            for i in range(N_kP) :
                write_string  = str(k_acc[i]) 
                for j in range(N_bands):  
                    pBand[k][j][i] -= E_Fermi
                    write_string += ','+ str(pBand[k][j][i]); 
                #    if Band[j][i]> Band_max: Band_max = pBand[0][j][i]
                #    if Band[j][i]< Band_min: Band_min = pBand[0][j][i] 
                write_string +='\n' 
                file0.write(write_string)
                
            file0.close()            
#------------------------------------------------------------------------------------------    
def INCAR_html_parameter():           
    global N_INCAR
    N_INCAR =0 
    name = []    
    file0 = open("./templates/INCAR_parameters.txt", 'r')
    while True:
        line = file0.readline()                   
        if not line : break
        else : line_list = line.split('=') ; name.append(line_list[0])            
        N_INCAR += 1
    file0.close()
    return name
#------------------------------------------------------------------------------------------        
def INCAR_html_explanation():      
    explanation = []
    file0 = open("./templates/INCAR_parameters.txt", 'r')
    while True:
        line = file0.readline()           
        if not line : break
        else : line_list = line.split('=') ; explanation.append(line_list[1])
    file0.close()
    return explanation    
#-----------------------------------------------------------------    
def INCAR_bullet():     
    global N_INCAR 
    row = [0]*N_INCAR

    bullet_is = [0,4,11,16,28,33,36,45,49]

    for i in range(9):
        row[bullet_is[i]] = 1
   
    return row   

    #-----------------------------------------------------------------    
def Check_Is_iter(string_check):  
    # here jemove only after of !   
    Is_iter =0
    CAR_i   =0
    CAR_iter_i =-1
    astring =-1
    
    for astring in string_check:
        if astring.find('{') !=-1 and astring.find('}') !=-1 : 

            Is_iter   += 1
            CAR_iter_i = CAR_i

            CAR_text_line = []
            
            CAR_text_line = astring.split('\n'); i=0
            for aline in CAR_text_line: CAR_text_line[i] = aline.split('!')[0]; i+=1   
            astring = "\n".join(CAR_text_line)     

         
        CAR_i += 1                           

    if Is_iter ==0: return [False, CAR_iter_i, astring] 
    else :          return [True , CAR_iter_i, astring] 

def Replace_iter(astring, replace_string):

    i_start = astring.find('{')
    i_end   = astring.find('}')    

    sel_string = astring[i_start+1: i_end]

    sel_string_list = sel_string.split()     
    sel_string = " ".join(sel_string_list)
    N_b = sel_string.count(' ') + 1
    
    for i in range(N_b):
        sel_string_list[i] = replace_string

    replace_string = " ".join(sel_string_list)

   
    astring = astring[0: i_start] + str(replace_string) + astring[i_end+1: ]

    return astring 
#------------------------------------------------------------------------------------------        
def atom_info_html_table():      
    #pt_atom_num = [[0]*18]*10
    pt_atom_num = [[0 for _ in range(18)] for _ in range(10)]         

    atom_p_pt = [ 
    [0,0],                                                                                                   [0,17],
    [1,0],[1,1],                                                                [1,12],[1,13],[1,14],[1,15],[1,16], [1,17],
    [2,0],[2,1],                                                                [2,12],[2,13],[2,14],[2,15],[2,16], [2,17],
    [3,0],[3,1],[3,2], [3,3],[3,4],[3,5],[3,6],[3,7],[3,8],[3,9],[3,10],[3,11], [3,12],[3,13],[3,14],[3,15],[3,16], [3,17],
    [4,0],[4,1],[4,2], [4,3],[4,4],[4,5],[4,6],[4,7],[4,8],[4,9],[4,10],[4,11], [4,12],[4,13],[4,14],[4,15],[4,16], [4,17],
    [5,0],[5,1],[5,2], 
                       [8,3],[8,4],[8,5],[8,6],[8,7],[8,8],[8,9],[8,10],[8,11],[8,12],[8,13],[8,14],[8,15],[8,16],
                       [5,3],[5,4],[5,5],[5,6],[5,7],[5,8],[5,9],[5,10],[5,11],[5,12],[5,13],[5,14],[5,15],[5,16],[5,17],
    [6,0],[6,1],[6,2], 
                       [9,3],[9,4],[9,5],[9,6],[9,7],[9,8],[9,9],[9,10],[9,11],[9,12],[9,13],[9,14],[9,15],[9,16],
                       [6,3],[6,4],[6,5],[6,6],[6,7],[6,8],[6,9],[6,10],[6,11],[6,12],[6,13],[6,14],[6,15],[6,16],[6,17]
    ]
   
    c_i = 1
    for atom_p in atom_p_pt:
        pt_atom_num[atom_p[0]][atom_p[1]] = c_i
        
        c_i += 1
    
    Atom_Name = ["",
		       "H" ,                                                                                 "He",
		       "Li","Be",                                                   "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
		       "Na","Mg",                                                   "Al","Si","P" ,"S" ,"Cl","Ar",
		       "K" ,"Ca","Sc", "Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
			   "Rb","Sr","Y" , "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
			   "Cs","Ba","La", 
                               "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
			                   "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
		       "Fr","Ra","Ac", 
                               "Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
			                   "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"  ]
    pseudo_type = [  ['LDA','PBE'], 
    ['LDA','PBE'],                                                                                                                                                                                                                                 ['LDA','PBE'],
    ['LDA','PBE'],['LDA','PBE'],                                                                                                                                             ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
    ['LDA','PBE'],['LDA','PBE'],                                                                                                                                             ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
    ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
    ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
    ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], 
                                              ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],
                                              ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
    ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], 
                                              ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],
                                              ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE']
    ]                               

    return pt_atom_num, Atom_Name, pseudo_type

def isNumber(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

#-----------------------------------------------------------------
basedir = os.path.abspath(os.path.dirname(__file__))
app = Flask(__name__)
app.config.update(
    UPLOADED_PATH=os.path.join(basedir, 'uploads'),
    # Flask-Dropzone config:
    #DROPZONE_ALLOWED_FILE_TYPE='image',
    #DROPZONE_ALLOWED_FILE_TYPE='.pdf',    
    DROPZONE_MAX_FILE_SIZE=1024,    
    DROPZONE_MAX_FILES=30,
    DROPZONE_TIMEOUT=5 * 60 * 1000  # set upload timeout to a large number, here is 5 minutes
)

app.config['DROPZONE_ALLOWED_FILE_CUSTOM'] = True
#app.config['DROPZONE_ALLOWED_FILE_TYPE'] = '.zip, .txt'
app.config['DROPZONE_ALLOWED_FILE_TYPE'] = '.zip, .txt'

dropzone = Dropzone(app)

@app.route("/", methods=['POST', 'GET'])
def Login():
    return render_template('Login_DB.html', F_initial = '0', User_text="", pass_text = "") 


@app.route("/DB", methods=['POST', 'GET'])
def start():

    print('start ----------------- ')
    
    result        = request.form
    curr_User     = result['User_info']
    curr_password = result['password_info']
    LoginID       = result['LoginID_info']

    if curr_User == "" or curr_password == "":        
        return render_template("Login_DB.html", F_initial = '1', User_text=curr_User, LoginID = LoginID)

    client = MongoClient('localhost', 27017)
    db = client['NewDBmat']

    User_DBdata = db.client['NewDBmat']["RUser"].find({'User': curr_User})
    
    F_exit = 0

    for di in User_DBdata: 
        F_exit = F_exit + 1        
        break
        
    if F_exit == 0: return render_template("Login_DB.html", F_initial = '1', User_text="\"" +curr_User+ "\"", LoginID = LoginID)

    Saved_password = di['User_password']    
    

    F_match = bcrypt.checkpw(curr_password.encode('utf-8'), Saved_password.encode('utf-8') )

    if F_match == True:    

        Dir_list = get_UserDir(curr_User) 

        User_Path      = di['User_Path']

        User_Paths     = []

        for (path, dirs, files) in os.walk(User_Path, followlinks=False):           

            User_Paths.append(path)

        User_Paths_txt = str(User_Paths)
        User_Paths_txt = User_Paths_txt.lstrip("[")
        User_Paths_txt = User_Paths_txt.rstrip("]")
        User_Paths_txt = User_Paths_txt.replace("'","")

        Key_List, User_List, Year_List, Journal_List, Model_List, program_List= get_DBinfo(curr_User) # The reason that 5 keywords are shown from start
        return render_template('NewLocalDB_ini.html', keys = Key_List, User_List = User_List, Year_List = Year_List, Journal_List = Journal_List, Model_List = Model_List, program_List = program_List, F_mod=0, N_db = 0, db_result = 0, curr_User = curr_User, User_Paths = User_Paths_txt, time_string = time_string, Dir_list= Dir_list) 
    else              : 
        return render_template('Login_DB.html', F_initial = '1', User_text="\"" +curr_User+ "\"", LoginID = LoginID)    

    

@app.route("/Update_DBWeb", methods=['POST', 'GET'])
def Update_DBWeb():

    result = request.form

    print('Update_DBWeb')
    print(result)

    #result = request.form
    curr_User = result['curr_User']    
    Year      = result['Year']
    Journal   = result['Journal']
    Model     = result['Model']
    Comment   = result['Comment']
    argout    = result['argout']

    selpath   = result['selpath']
    selpath   = selpath.strip()

    F_mod = 0

    

    update_DB(curr_User,Year,Journal,Model,Comment,argout,selpath)

    Key_List, User_List, Year_List, Journal_List, Model_List, program_List = get_DBinfo(curr_User)       

    Dir_list = get_UserDir(curr_User)

    F_mod = 1 

    DB_updateWeb_html = "temp_DB_updateWeb"+ time_string +".html"

    return render_template(DB_updateWeb_html, keys = Key_List, User_List = User_List, Year_List = Year_List, Journal_List = Journal_List, Model_List = Model_List, program_List = program_List, F_mod = F_mod, selpath = selpath, Dir_list = Dir_list) 

#-----------------------------------------------------------------
@app.route("/NewLocalDB", methods=['POST', 'GET'])
def NewLocalDB():
        
    result = request.form
    
    static_path= "/home3/DB_devel/static/"
    html_path  = "htmls/"
    static_html_path = static_path + html_path
    
    if str(result).find("'search_text'") != -1:    
               
        client = MongoClient('localhost', 27017)

        db = client['NewDBmat']
        
        keywords =['User','Journal', 'Program', 'Model', 'Calculation_type', 'Directory_level']
        
        key_val = {}

        for ki in keywords:                   
            val = result.getlist(ki)[0]
            
            if val !='': key_val[ki] = val                
        
        search_text = result['search_text']
        
        if search_text !='':
            #for ki in keywords:

            #    key_search_0 = {ki:{"$regex":"EGF"},"Info_path":{"$regex":"EGF"} }

            key_search_0 = {"$regex":search_text }            
            key_val["Info_path"] = key_search_0
            
        DFT_DBdata = db.client["DBmat"]["Group"].find(key_val)

        N_db = db.client["DBmat"]["Group"].count(key_val)

        # --- create --- iframe -------------------------------
        center_html = make_html.center_string(DFT_DBdata, search_text)

        center_local_url = html_path + "temp_center_"+str(int(time.time()*10000000))+".html"

        center_url = static_path + center_local_url                
        command = "rm -f " + static_html_path + "temp_center*.html"; os.system(command)

        p_file = open(center_url, 'w')
        p_file.write(center_html)
        p_file.close()


        hidden_info_field = result['hidden_info_field']
       

        N_field   = []
        selected_i= []

        field_key = []

                

        field_info0 = hidden_info_field.split('/')

      

        for i, f0_i in enumerate(field_info0):
            if i==0:continue
            field_info1 = f0_i.split('@@@')
        
            N_field.append(field_info1[1])
            selected_i.append(field_info1[2])

            field_info2 = field_info1[3].split('_@@')

            a_field = []

            for j, f1_i in enumerate(field_info2):                
                if j==0:continue
                a_field.append(f1_i)

            field_key.append(a_field)

        #---------------------------------------------------------------------

        return render_template("NewLocalDB.html", center_url = center_local_url, DB_key=result, N_db=N_db, N_field=N_field, selected_i=selected_i, field_key=field_key )

    
   
#-----------------------------------------------------------------
@app.route("/SearchDB", methods=['POST', 'GET'])
def SearchDB():
    result = request.form

    client = MongoClient('localhost', 27017)

    db = client['NewDBmat']

    field_info = result['field_info']
       
    arr_field_info = field_info.split('/')

    key_val = {}

    for i, field_i in enumerate(arr_field_info):
        if i==0:continue
        field_info1 = field_i.split('@@@')            
        key_val[field_info1[1]] = field_info1[2]
    
    
    
    search_text = result['search_text']
    

    if search_text != "": 

        arr_search_text = search_text.split(" ")

        if len(arr_search_text) > 1:

            temp_list = []

            for atext in arr_search_text:            
                temp_list.append({'$regex':atext})
            
            key_search_0 = tuple(temp_list)

        else :

            key_search_0 = {'$regex':arr_search_text[0]}


        
        key_val['formula'] = key_search_0

        print('key_val', key_val)
        
    #key_val['formula'] = "C60"

    db_list = db.client["NewDBmat"]["DFT"].find(key_val)

    print('db_list', db_list)
    N_db    = db.client["NewDBmat"]["DFT"].count_documents(key_val)

    print('N_db', N_db)

    
    db_result = ""

    for i, di in enumerate(db_list):	

        #db_result.append(di)
        str_di = str(di)
        str_di = str_di.lstrip("{")
        str_di = str_di.rstrip("}")

        str_di = str_di.replace("\'","")

        db_result += str_di + "@@@"
    
    #str_db_result = "\"" + str(db_result) + "\""
    #str_db_result = str(db_result)
    #str_db_result = str_db_result.lstrip("[")
    #str_db_result = str_db_result.rstrip("]")

    #str_db_result = str_db_result.replace("\'","\"")
    
    

    DB_search_html = "temp_DB_search"+ time_string +".html"
        
    return render_template(DB_search_html, F_mod=1, N_db = N_db, db_result = db_result)

@app.route('/uploads', methods=['POST', 'GET'])
def uploads():
    result = request.form

    User = result['User_info']
    
    if request.method != 'POST': return 0

    f = request.files.get('file')
    file_path = os.path.join(app.config['UPLOADED_PATH'], f.filename)
    f.save(file_path)
    

    command = "rm -rf ./uploads/" + User;    os.system(command)  # remove old data
    command = "unzip -qq " + file_path + " -d /home3/DB_devel/uploads/" + User; os.system(command)    
    command = "rm -f " + file_path;     os.system(command)


    Dropzone_html = "temp_Dropzone"+ time_string +".html"

    return render_template(Dropzone_html)
    

#-----------------------------------------------------------------
@app.route("/Logout_DB", methods=['POST', 'GET'])
def Logout():
    result = request.form   

    #User = result.getlist('User_info')[0] 

    return render_template('Login_DB.html')

#-----------------------------------------------------------------
@app.route("/Register_DB", methods=['POST', 'GET'])
def Register():
    
    return render_template("Register_DB.html", Message_return = "", User_text="", password_text="")

#-----------------------------------------------------------------
@app.route("/Register_ID_DB", methods=['POST', 'GET'])
def Register_ID_DB():
    result = request.form
    
    User = result['User_info']
    User_password = result['password_info']
    User_Path = result['Path_info']
      
    client = MongoClient('localhost', 27017)
    db = client['NewDBmat']
    User_DBdata = db.client['NewDBmat']["RUser"].find({'User': User})

    F_exit = 0

    for di in User_DBdata: F_exit = F_exit + 1

    if F_exit > 0:        
        return render_template("Register_DB.html", Message_return = "ID that already exists", Message_color = "red", User_text=User, password_text=User_password, Path_text=User_Path)

    
    passwordSalt = bcrypt.gensalt()

    User_password = User_password.encode('utf-8')

    passwordHash = bcrypt.hashpw(User_password, passwordSalt)
    
    insertPasswordHash = passwordHash.decode() 
    
    db.client["NewDBmat"]["RUser"].insert_one({ 'User' : User,  'User_password' : insertPasswordHash,  'LastSN' : 0,  'User_Path' : User_Path })

    work_dir  = "./static/DB_Group/" + User
    
    command_string = "mkdir " + work_dir ;  os.system(command_string)

    Download_dir  = "./static/Downloads/" + User
    
    command_string = "mkdir " + Download_dir ;  os.system(command_string)

    return render_template("Register_DB.html", Message_return = 'User was successfully created', Message_color = "blue")

#-----------------------------------------------------------------
@app.route("/Download", methods=['POST', 'GET'])
def Download():

    result = request.form

    New_path_file = result['NewPath_info']

    curr_User = New_path_file.split("__")[0]
    
    target_dir = "./static/DB_Group"

    targetUser_dir = target_dir + "/" + curr_User
    

    New_path = targetUser_dir + "/" + New_path_file

    print(New_path)

    Download_path = "/home3/DB_devel/static/Downloads/" + curr_User

    zip_path = Download_path + '/' + New_path_file + ".zip"                

    # remove all earlier zip before making new zip  -------------------------                
    command = "rm -f " + Download_path + "/*"; os.system(command)

    # make zip download -------------------------                
    command = "cd " + New_path + ";zip " + zip_path + " ./*"; os.system(command)


    Download_html = "temp_DB_download"+ time_string +".html"

    return render_template(Download_html, F_mod = 1)

    

def update_DB(curr_User,Year,Journal,Model,Comment,argout,selpath):
    print("Update DB----------------")

    if selpath =="":        
        work_dir = "./uploads"
        workUser_dir = work_dir + "/" + curr_User
    else :                
        workUser_dir = selpath

    
    target_dir = "./static/DB_Group"
    targetUser_dir = target_dir + "/" + curr_User

    client = MongoClient('localhost', 27017)
    db = client['NewDBmat']
    User_DBdata = db.client['NewDBmat']["RUser"].find({'User': curr_User})

    for di in User_DBdata:        
        Last_SN = di['LastSN']
        break

    
    SN = int(Last_SN) + 1

    for (path, dirs, files) in os.walk(workUser_dir, followlinks=False):    

        #print('path ', path)       

        is_DATA   = False

        is_DFT = {'VASP': False, 'SIESTA': False, 'TranSIESTA': False, 'MSDFT': False, 'Transport': False}        

        program   = ""

        
        fdflog_dir_file   = path + '/fdf-*.log' 
        fdflog_file_list = glob(fdflog_dir_file)
        
        if fdflog_file_list != []:
            is_DATA            = True
            fdf_log, program   = get_fdf_log(fdflog_file_list[0])
            is_DFT[program]    = True

        OUTCAR_file     = path + "/OUTCAR"
        Is_OUTCAR_file  = os.path.isfile(OUTCAR_file)

        if Is_OUTCAR_file == True:            
            is_DATA         = True
            program         = "VASP"
            is_DFT[program] = True


        if is_DATA == False: continue               
        

        # Data dir ------------------------------------------------------------------
        # common 

        arr_path = path.split("/")
             
        if len(arr_path) > 3: 
            Comment_dir = '_'.join(arr_path[4:])            
        else:
            Comment_dir = "No Comment_dir"                

        str_SN = str(SN)
        str_SN = str_SN.zfill(5)

        SN += 1

        # copy data to target -------------------------
        New_path_file = curr_User + "__" + str_SN + "__" + Year + "__" + Journal + "__" + Model + "__" + program + "__" + Comment
        New_path_dir  = targetUser_dir + "/" + program
        New_path      = New_path_dir + "/" + New_path_file
        

        Is_dir  = os.path.isdir(New_path_dir)

        if Is_dir == False:
            command = "mkdir " + New_path_dir; os.system(command)

        #New_path = targetUser_dir + "/" + New_path_file

        command = "cp -r " + path + " " + New_path; os.system(command)

        # SIESTA ------------------------------------------------------------------------
        if program == "SIESTA" or program == "TranSIESTA" or program == "MSDFT" or program == "Transport":       

            #  ---- Generate visualizing
            command = '/home3/DB_devel/Utils/post_db ' + New_path   ; os.system(command)

            merge_fdf = get_merge_fdf(New_path)

            argout = argout.strip()
  
            addDB_SIESTA(New_path,db, files,Comment_dir, merge_fdf, fdf_log, argout)           
        
        # VASP ------------------------------------------------------------------------
        if program == "VASP":
            
            command = '/home3/DB_devel/Utils/post_VASP ' + New_path   ; os.system(command)            
            addDB_VASP(New_path,db, files,Comment_dir, argout)                  
            #Read_xml(New_path,1)             
            #Make_graph_png(New_path)       
    
    Last_SN = SN - 1
    db.client['NewDBmat']["RUser"].update_one({'User': curr_User}, {"$set":{"LastSN":Last_SN}})

def Make_graph_png(work_dir):
    global Is_spin_dn

    plt.clf()         
    filename_string  = work_dir + "/tDOS.dat" 
    graph_data = pd.read_csv(filename_string)
    
    graph_x = graph_data['energy']; graph_y = graph_data['total_up']; plt.plot(graph_x, graph_y)

    
    if Is_spin_dn: graph_y = -graph_data['total_dn']; plt.plot(graph_x, graph_y)
    
    if Is_pDOS: 
        
        filename_string  = work_dir + "/pDOS.dat" 
        graph_data = pd.read_csv(filename_string)
        

        graph_x = graph_data[DOS_title[0]]
        
        for i in range(1, N_p_title):                            
            graph_y =  graph_data[DOS_title[i]+'_up']; plt.plot(graph_x, graph_y)               
            if Is_spin_dn :
                graph_y = -graph_data[DOS_title[i]+'_dn']; plt.plot(graph_x, graph_y)



    #plt.xlim([-10, 10])
    

    plt.xlabel('E(eV)')
    plt.ylabel('DOS')
    plt.legend(loc='upper left')
    plt.title("Density of state")

    graph_process(work_dir)     

def get_current(New_path, argout):   

    arg_file = New_path + "/" + argout     
    Is_arg_file  = os.path.isfile(arg_file)

    if Is_arg_file == False:
        return "No", "No", "No", "No"

    p_arg_file = open(arg_file,'r')
    contents = p_arg_file.read()
    p_arg_file.close()

    arr_contents = contents.split("\n") 

    On_Currents = False

    Currents_i = 0

    Current      = ""
    Current_unit = ""
    Power        = ""
    Power_unit   = ""

    for a_line in arr_contents:
        a_line = a_line.strip()
        arr_a_line = a_line.split(" ")        

        if On_Currents == True and Currents_i ==2:

            arr_a_line = a_line.split("/")           

            Power_info = arr_a_line[2].strip()

            arr_Power_info = Power_info.split(" ")
            Power      = arr_Power_info[0]
            Power_unit = arr_Power_info[1]
            break

        if On_Currents == True and Currents_i ==1:

            arr_a_line = a_line.split("/")            

            Current_info = arr_a_line[2].strip()

            arr_Current_info = Current_info.split(" ")
            Current      = arr_Current_info[0]
            Current_unit = arr_Current_info[1]

            Currents_i = 2


        if arr_a_line[0] == "Currents": 
            On_Currents = True
            Currents_i = 1

    
    return Current, Current_unit, Power, Power_unit


def get_fdf_log(fdflog_file):   

    p_fdflog_file = open(fdflog_file,'r')
    contents = p_fdflog_file.read()
    p_fdflog_file.close()

    arr_contents = contents.split("\n")

    New_contents = ''

    program_type = "SIESTA"

    for a_line in arr_contents:
        i_sharp = a_line.find("#")
        if i_sharp == 0 : continue
        if i_sharp !=-1 : a_line = a_line[0:i_sharp] 

        a_line = ' '.join(a_line.split())        

        New_contents += a_line +"\n"    

        arr_a_line = a_line.split(" ")

        if arr_a_line[0] == "MSDFT"            : program_type = "MSDFT"
        if arr_a_line[0] == "TS.SolutionMethod": program_type = "TranSIESTA"
        if arr_a_line[0] == "TS.Voltage"       : program_type = "Transport"
        if arr_a_line[0] == "MSCDFTBias"       : program_type = "MSDFT"

    
    return New_contents, program_type

def get_merge_fdf(New_path):
    fdf_dir       = New_path
    fdf_file      = fdf_dir + '/*.fdf' 
    fdf_file_list = glob(fdf_file)

    # Fine the main fdf file ------------------------------------------------
    Main_contents = ''
    Sub_contents  = {}
    # Arrange files ------------------------------------------------
    for a_fdf in fdf_file_list:
        file0 = open(a_fdf, 'r')
        contents = file0.read()
        file0.close()
        Is_include = False
        
        arr_contents = contents.split('\n')    
        New_contents = ''

        for a_line in arr_contents:
            
            i_sharp = a_line.find("#")
            if i_sharp == 0 : continue
            if i_sharp !=-1 : a_line = a_line[0:i_sharp] 

            a_line = ' '.join(a_line.split())        
            if a_line=='' : continue

            arr_a_line = a_line.split(" ")

            if arr_a_line[0] == '%include':
                Is_include = True
                Main_fdf = a_fdf
        
            New_contents += a_line + "\n"
        
        #New_contents = New_contents.rstrip("\n")

        
        if Is_include == True or len(fdf_file_list) == 1 :             
            Main_fdf            = a_fdf
            Main_contents       = New_contents        
        else                  : 
            Sub_contents[a_fdf] = New_contents

    #--------------------------------------------------------

    Main_fdf = fdf_dir + Main_fdf 

    arr_contents = Main_contents.split('\n')

    New_contents = ''

    for a_line in arr_contents:
        arr_a_line   = a_line.split(" ")
        if arr_a_line[0] == '%include':              
            a_file = fdf_dir + "/"+ arr_a_line[1]
            a_line = Sub_contents[a_file]       

        New_contents += a_line + '\n'    

    return New_contents

def addDB_VASP(path,db,files,Comment_dir, argout):
    path_dir_arr = path.split('/')
    str_dir = path_dir_arr[-1]
    arr_str_dir = str_dir.split("__")            

    DB_item = {}
    DB_item['User']     = arr_str_dir[0]
    DB_item['DataSN']   = arr_str_dir[1]

    DB_item['Year']     = arr_str_dir[2]
    DB_item['Journal']  = arr_str_dir[3]
    DB_item['Model']    = arr_str_dir[4]

    DB_item['program']  = arr_str_dir[5]        
    DB_item['Comment']  = arr_str_dir[6]
    DB_item['Comment_dir']  = Comment_dir
    DB_item['Arg_out']  = argout

    vasprun_file = path + "/vasprun.xml"

    vasp = vasprun(vasprun_file)

    #'/home/joonho/research/H2-Pt-C60-x/C60/1FC60Pt66H2/vasprun.xml'
    #/home/joonho/research/H2-Pt-C60-x/C60/1FC60Pt66H2


    for key in vasp.values.keys():    

        

        key = key.strip()

        if key =="name_array" : continue
        if key =="composition": continue
        if key =="elements"   : continue
        if key =="finalpos"   : continue

        if key =="kpoints"    : continue


        if str(type(vasp.values[key])) == "<class 'dict'>": 
            temp_pseudo_string = ""
            F_pseudo  =0
            for akey in vasp.values[key].keys():

                akey = akey.strip()

                if key =="calculation" :
                    if akey =="tdos": continue   
                    if akey =="eband_eigenvalues": continue      
                    if akey =="force": continue
                    if akey =="pdos": continue
                    if akey =="projected": continue
                    if akey =="born_charges": continue
                    if akey =="hessian": continue
                    if akey =="normal_modes_eigenvalues": continue
                    if akey =="normal_modes_eigenvectors": continue
                    if akey =="epsilon_ion": continue
                    if akey =="pion": continue
                    if akey =="psp1": continue
                    if akey =="psp2": continue
                    if akey =="pelc": continue
                    if akey =="force": continue

                if str(type(vasp.values[key][akey])) == "<class 'dict'>": 
                    for aakey in vasp.values[key][akey].keys():    
                        aakey = aakey.strip()                
                        
                        if str(type(vasp.values[key][akey][aakey])) == "<class 'dict'>": 
                            for aaakey in vasp.values[key][akey][aakey].keys():                    
                                aaakey = aaakey.strip()
                                
                                DB_item[aaakey] = vasp.values[key][akey][aakey][aaakey]

                        elif str(type(vasp.values[key][akey][aakey])) == "<class 'list'>":                

                            if   np.array(vasp.values[key][akey][aakey]).ndim == 1:                         
                                temp_string = remove_side(vasp.values[key][akey][aakey], "[", "]" )       
                                

                                DB_item[aakey] = temp_string
                            
                            elif np.array(vasp.values[key][akey][aakey]).ndim == 2:                    

                                temp_string = ""
                                for item in vasp.values[key][akey][aakey]:                                    
                                    temp_string += remove_side(item, "[", "]" ) + "<br>"
                                temp_string = temp_string.rstrip("<br>")                                
                                
                                temp_string = "BlockInfo*@@*" + temp_string.replace(","," ")
                                
                                DB_item[aakey] = temp_string

                        else:
                            temp_string = str(vasp.values[key][akey][aakey])
                            
                            temp_string = temp_string.replace(","," ")
                            
                            DB_item[aakey] = temp_string

                elif str(type(vasp.values[key][akey])) == "<class 'list'>":                

                    
                    if key == 'pseudo_potential': 
                        F_pseudo += 1                    
                        temp_pseudo_string += remove_side(vasp.values[key][akey], "[", "]" ) + "<br>"
                        if F_pseudo == 3 : 
                            
                            temp_pseudo_string  = temp_pseudo_string.rstrip("<br>")
                            temp_pseudo_string  = "BlockInfo*@@*" + temp_pseudo_string.replace(","," ")     
                            

                            DB_item[key] = temp_pseudo_string
                            F_pseudo = 0
                        
                    else :    

                        if   np.array(vasp.values[key][akey]).ndim == 1:

                            temp_string = remove_side(vasp.values[key][akey], "[", "]" )
                            
                            
                            temp_string = temp_string.replace(","," ")
                            
                            DB_item[akey] = temp_string

                        elif np.array(vasp.values[key][akey]).ndim == 2:                              

                            temp_string = ""
                            for item in vasp.values[key][akey]:                                
                                temp_string += remove_side(item, "[", "]" ) + "<br>"
                            temp_string = temp_string.rstrip("<br>")                            
                            temp_string = "BlockInfo*@@*" + temp_string.replace(","," ")                                                        
                            DB_item[akey] = temp_string

                else:

                    temp_string = str(vasp.values[key][akey])

                    i_Exclamation = temp_string.find("!")
                    
                    if i_Exclamation !=-1 : temp_string = temp_string[0:i_Exclamation] 

                    temp_string = temp_string.strip()                                         
                    temp_string = temp_string.replace(","," ")
                    
                    DB_item[akey] = temp_string

        elif str(type(vasp.values[key])) == "<class 'list'>":

            if np.array(vasp.values[key]).ndim == 1:
                temp_string = remove_side(vasp.values[key], "[", "]" )                
                temp_string = temp_string.replace(","," ")                
                DB_item[key] = temp_string

            elif np.array(vasp.values[key]).ndim == 2:            

                temp_string = ""
                for item in vasp.values[key]:                    
                    temp_string += remove_side(item, "[", "]") + "<br>"
                temp_string = temp_string.rstrip("<br>")                
                temp_string = "BlockInfo*@@*" + temp_string.replace(","," ")                                
                DB_item[key] = temp_string
        
        else:

            temp_string = str(vasp.values[key])
            i_Exclamation = temp_string.find("!")            
            if i_Exclamation !=-1 : temp_string = temp_string[0:i_Exclamation] 
            temp_string = temp_string.strip()                   
            temp_string = temp_string.replace(","," ")                   
                
            DB_item[key] = temp_string
            
    db.client["NewDBmat"]["DFT"].insert_one(DB_item)

    N_files = len(files)
    
    File_list =''
    for file_name in files:               

        path_filename = path + "/" + file_name
        
        if os.path.islink(path_filename) == True: continue

        File_list += file_name + "\n"


    DirInfoFile = path + "/DB_info.inf"

    p_file = open(DirInfoFile, 'w')

    write_string="@File_list " + str(N_files) +"\n" + File_list     
    p_file.write(write_string)    
    p_file.close()

def addDB_SIESTA(path,db,files,Comment_dir, merge_fdf, fdf_log, argout):

    path_dir_arr = path.split('/')

    str_dir = path_dir_arr[-1]

    arr_str_dir = str_dir.split("__")            

    DB_item = {}
    DB_item['User']     = arr_str_dir[0]
    DB_item['DataSN']   = arr_str_dir[1]

    DB_item['Year']     = arr_str_dir[2]
    DB_item['Journal']  = arr_str_dir[3]
    DB_item['Model']    = arr_str_dir[4]

    DB_item['program']  = arr_str_dir[5]        
    DB_item['Comment']  = arr_str_dir[6]
    DB_item['Comment_dir']  = Comment_dir
    DB_item['Arg_out']  = argout

    if DB_item['program'] == "Transport":
        Current, Current_unit, Power, Power_unit = get_current(path, argout)
        DB_item['Current']      = Current
        DB_item['Current_unit'] = Current_unit
        DB_item['Power']        = Power
        DB_item['Power_unit']   = Power_unit

    # merge fdf ---------------------------------------------------

    block_ON = False        

    arr_merge_fdf = merge_fdf.split("\n")

    for fdf_line in arr_merge_fdf:        

        if fdf_line == "": continue
        
        arr_fdf_line = fdf_line.split(" ")
        Field_name = arr_fdf_line[0].lower()        

        if Field_name == "%endblock" : 
            block_ON = False
            BlockInfo = BlockInfo.rstrip("<br>")
            DB_item[block_type] = BlockInfo
            continue

        if block_ON == True: 
            BlockInfo += fdf_line + "<br>"
            continue

        if Field_name == "%block": 
            block_type = arr_fdf_line[1].lower()
            block_type = block_type.replace(".","_")
            block_type = block_type.replace("$","_")

            block_ON = True
            BlockInfo = "BlockInfo*@@*"            
            continue
        
    
    # fdf log ---------------------------------------------------


    arr_fdf_log = fdf_log.split("\n")

    for fdf_line in arr_fdf_log:        

        if fdf_line == "": continue

        i_sharp = fdf_line.find("#")
        if i_sharp == 0 : continue
        if i_sharp !=-1 : fdf_line = fdf_line[0:i_sharp] 

        fdf_line = ' '.join(fdf_line.split())                
        arr_fdf_line = fdf_line.split(" ")

        if len(arr_fdf_line) == 1 : continue       


        Field_name = arr_fdf_line[0]
        Field_name = Field_name.replace(".","_")
        Field_name = Field_name.replace("$","_")

        if Field_name == "%endblock"   : continue
        if isNumber(Field_name) == True: continue
        if Field_name == "%block"      : continue
        
        if Field_name == 'SystemName':        
            DB_item[Field_name] = ' '.join(arr_fdf_line[1:])
            continue
        
        
        DB_item[Field_name] = arr_fdf_line[1]
        
        
        if len(arr_fdf_line) == 3 :        

            temp_Field_name = Field_name + "_unit"

            DB_item[temp_Field_name] = arr_fdf_line[2]    
            

    # OUTVARS.yml  ---------------------------------------------------

    OUTVARS_file = path + '/OUTVARS.yml'
    Is_OUTVARS  = os.path.isfile(OUTVARS_file)
    if Is_OUTVARS == True:
        OUTVARS_energies_ON = False

        p_OUTVARS_file = open(OUTVARS_file, 'r')  
        while True:
            read_string = p_OUTVARS_file.readline()
            if not read_string : break                
            read_string = read_string.strip()
            read_string = ' '.join(read_string.split())
            arr_read_string = read_string.split(" ")

            if OUTVARS_energies_ON == True:
                if arr_read_string== [''] : 
                    OUTVARS_energies_ON = False
                    break

                
                E_keyword = arr_read_string[0].rstrip(":")
                DB_item[E_keyword] = arr_read_string[1].strip()
                
                

            if arr_read_string[0] =="energies:":
                OUTVARS_energies_ON = True

            if arr_read_string[0] =="version:":
                Siesta_version = arr_read_string[1].strip("\"")
                DB_item['version'] = Siesta_version.strip()
                
        
        p_OUTVARS_file.close()
    else:
        print("There is no OUTVARS.yml")

    # BASIS_ENTHALPY  ---------------------------------------------------

    BASIS_ENTHALPY_file = path + '/BASIS_ENTHALPY'
    Is_BASIS_ENTHALPY_file  = os.path.isfile(BASIS_ENTHALPY_file)
    if Is_BASIS_ENTHALPY_file == True:

        p_BASIS_ENTHALPY_file = open(BASIS_ENTHALPY_file, 'r')  
        read_string = p_BASIS_ENTHALPY_file.readline()
        read_string = read_string.strip()

        DB_item['BASIS_ENTHALPY'] = read_string                
                        
        p_BASIS_ENTHALPY_file.close()
    else:
        print("There is no BASIS_ENTHALPY")

    # BASIS_HARRIS_ENTHALPY  ---------------------------------------------------

    BASIS_HARRIS_ENTHALPY_file = path + '/BASIS_HARRIS_ENTHALPY'
    Is_BASIS_HARRIS_ENTHALPY_file  = os.path.isfile(BASIS_HARRIS_ENTHALPY_file)
    if Is_BASIS_HARRIS_ENTHALPY_file == True:

        p_BASIS_HARRIS_ENTHALPY_file = open(BASIS_HARRIS_ENTHALPY_file, 'r')  
        read_string = p_BASIS_HARRIS_ENTHALPY_file.readline()
        read_string = read_string.strip()

        DB_item['BASIS_HARRIS_ENTHALPY'] = read_string                
                        
        p_BASIS_HARRIS_ENTHALPY_file.close()
    else:
        print("There is no HARRIS_ENTHALPY")
        
    # FORCE_STRESS  ---------------------------------------------------

    FORCE_STRESS_file = path + '/FORCE_STRESS'
    Is_FORCE_STRESS_file  = os.path.isfile(FORCE_STRESS_file)

    if Is_FORCE_STRESS_file == True:

        p_FORCE_STRESS_file = open(FORCE_STRESS_file, 'r')          
        read_string = p_FORCE_STRESS_file.readline()        
        read_string = read_string.strip()        
        DB_item['FORCE_STRESS'] = read_string                                        
        p_FORCE_STRESS_file.close()
    else:
        print("There is no FORCE_STRESS")

    # TIMES  ---------------------------------------------------

    TIMES_file = path + '/TIMES'
    Is_TIMES_file  = os.path.isfile(TIMES_file)
    if Is_TIMES_file == True:

        F_Nnode = False
        F_ETime = False

        p_TIMES_file = open(TIMES_file, 'r')  
        while True:
            read_string = p_TIMES_file.readline()
            if not read_string : break                

            read_string = read_string.strip()
            

            if "timer: Number of nodes =" in read_string:
                read_string = ' '.join(read_string.split())
                arr_read_string = read_string.split("=")
                DB_item['N_node'] = arr_read_string[1].strip()
                F_Nnode = True

            if "timer: Total elapsed wall-clock time (sec)" in read_string :
                read_string = ' '.join(read_string.split())
                arr_read_string = read_string.split("=")
                DB_item['Elapsed_time'] = arr_read_string[1].strip()
                F_ETime = True
                break
            

        if F_Nnode == False:    print("There is no timer: Number of nodes")
        if F_ETime == False:    print("There is no Total elapsed wall-clock time")

    else:
        print("There is no TIMES")


    db.client["NewDBmat"]["DFT"].insert_one(DB_item)

    N_files = len(files)
    
    File_list =''
    for file_name in files:
        

        path_filename = path + "/" + file_name
        
        if os.path.islink(path_filename) == True: continue

        File_list += file_name + "\n"


    DirInfoFile = path + "/DB_info.inf"

    p_file = open(DirInfoFile, 'w')

    write_string="@File_list " + str(N_files) +"\n" + File_list     
    p_file.write(write_string)    
    p_file.close()


def get_UserDir(curr_User):

    User_Path = "/home3/DB_devel/static/DB_Group/" + curr_User

    User_dirs = ""

    for (path, dirs, files) in os.walk(User_Path, followlinks=False):                  

        arr_path = path.split("/")

        print('arr_path')
        print(arr_path)

        if len(arr_path) == 7:           

            for dir in dirs:
                User_dirs += dir + ","
            #break    

    User_dirs = User_dirs.rstrip(",")

    print('User_dirs')
    print(User_dirs)

    return User_dirs



def get_DBinfo(curr_User):
    
    client = MongoClient('localhost', 27017)
    db     = client['NewDBmat']
    cursor = db.client["NewDBmat"]["DFT"].find({})

    exc_field = ["_id", "DataSN", "User", "Year", "Journal", "Model", "program"]

    keylist = []

    for document in cursor: 

        #keys = str(document.keys())
        keylist = keylist + list(document.keys())

        keylist = list(set(keylist))

    set_keylist   = set(keylist)
    set_exc_field = set(exc_field)

    keylist = list(set_keylist-set_exc_field)

    keylist.sort()

    Key_List = str(keylist) 

    Key_List = Key_List.lstrip("[")
    Key_List = Key_List.rstrip("]")
    Key_List = Key_List.replace("'","")    
    
    Key_List = "DataSN,User,Year,Journal,Model,program," + Key_List

    fixed_field = ["User", "Year", "Journal", "Model", "program"]

    member_List = []

    for a_field in fixed_field:    
        field_members   = db.client["NewDBmat"]["DFT"].find().distinct(a_field)
        field_members.sort()

        members_str = str(field_members) 

        members_str = members_str.lstrip("[")
        members_str = members_str.rstrip("]")
        members_str = members_str.replace("'","")           
        member_List.append(members_str)

    return Key_List, member_List[0], member_List[1], member_List[2], member_List[3], member_List[4]


if __name__ == "__main__":        

    command = "rm -f /home3/DB_devel/static/temp_*.html"; os.system(command)


    time_string = str(int(time.time()*1000000))

    

    NewLocalDB_ini_html = make_html.create_NewLocalDB_ini()

    DB_search           = make_html.create_DB_search(time_string)
    DB_updateWeb        = make_html.create_DB_updateWeb(time_string)
    Dropzone_html       = make_html.create_Dropzone_html(time_string)
    DB_download         = make_html.create_DB_download(time_string)

    Login_DB_html       = make_html.create_Login_DB_html()
    Logout_DB_html      = make_html.create_Logout_DB_html()
    Register_DB_html    = make_html.create_Register_DB_html()

    




    # for new.gnubands
    #command = "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/intel/Compiler/18.0.4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib"
    #os.system(command)

    

    


    app.run(host="0.0.0.0", port=80)

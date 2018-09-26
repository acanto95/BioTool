#!/usr/bin/env python
import sys,os
import re
import numpy as np
import urllib.request
import plotly.graph_objs as go
import plotly.offline as ply
import tkinter as tk
from tkinter import *


def show_text():


   nombre= entry_text.get()
   nombrearchivo = nombre + ".pdb"
   connection=urllib.request.urlopen("https://files.rcsb.org/view/"+nombre+".pdb")
   output=connection.read().decode('utf-8')


   # Write data to file
   filename = nombrearchivo
   file2 = open(filename, 'w+')
   file2.write(str(output))
   file2.close()

   aa3to1={
      'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
      'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
      'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
      'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
      'MSE':'M',
   }

   ca_pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
   for pdb_file in nombrearchivo:
       filename=os.path.basename(pdb_file).split('.')[0]
       chain_dict=dict()
       chain_list=[]

       fp=open(nombrearchivo,'rU')
       for line in fp.read().splitlines():
           if line.startswith("ENDMDL"):
               break
           match_list=ca_pattern.findall(line)
           if match_list:
               resn=match_list[0][0]+match_list[0][2]
               chain=match_list[0][1]+match_list[0][3]
               if chain in chain_dict:
                   chain_dict[chain]+=aa3to1[resn]
               else:
                   chain_dict[chain]=aa3to1[resn]
                   chain_list.append(chain)
       fp.close()


       data= ''
       for chain in chain_list:
          data = ('M'+chain + chain_dict[chain])



   Asp = float(data.count('D'))
   Thr = float(data.count('T'))
   Esr = float(data.count('S'))
   Glu = float(data.count('E'))
   Pro = float(data.count('P'))
   Gly = float(data.count('G'))
   Ala = float(data.count('A'))
   Cys = float(data.count('C'))
   Vla=  float(data.count('V'))
   Met = float(data.count('M'))
   Ile = float(data.count('I'))
   Leu = float(data.count('L'))
   Tyr = float(data.count('Y'))
   Phe = float(data.count('F'))
   His = float(data.count('H'))
   Lys = float(data.count('K'))
   Arg = float(data.count('R'))
   Trp = float(data.count('W'))
   Gln = float(data.count('Q'))
   Asn = float(data.count('N'))


   Asp = ((Asp*100)/424)
   Thr = ((Thr*100)/424)
   Esr = ((Esr*100)/424)
   Glu = ((Glu*100)/424)
   Pro = ((Pro*100)/424)
   Gly = ((Gly*100)/424)
   Ala = ((Ala*100)/424)
   Cys = ((Cys*100)/424)
   Vla = ((Vla*100)/424)
   Met = ((Met*100)/424)
   Ile = ((Ile*100)/424)
   Leu = ((Leu*100)/424)
   Tyr = ((Tyr*100)/424)
   Phe = ((Phe*100)/424)
   His = ((His*100)/424)
   Lys = ((Lys*100)/424)
   Arg = ((Arg*100)/424)
   Trp = ((Trp*100)/424)
   Gln = ((Gln*100)/424)
   Asn = ((Asn*100)/424)


   n = 201
   x = np.linspace(0, 2.0*np.pi, n)
   y1 = np.sin(x)
   y2 = np.cos(x)
   y3 = y1 + y2

   trace1 = go.Bar(
       x=['GLY(G): Glycine', 'ALA(A): Alanine', 'VAL(V): Valine','LEU(L): Leucine','ILE(I): Isoleucine','MET(M): Methionine','PHE(F): Phenylalanine','TRP(W): Tryptophan','PRO(P): Proline'],
       y=[Gly, Ala, Vla, Leu,Ile,Met,Phe,Trp,Pro],
       name='Group A: Nonpolar Amino Acids (Hydrophobic)'
   )
   trace2 = go.Bar(
       x=['SER(S): Serine', 'THR(T): Threonine', 'CYS(C): Cysteine','TYR(Y): Tyrosine','ASN(N): Asparagine','GLN(Q): Glutamine'],
       y=[Esr, Thr, Cys,Tyr, Asn, Gln],
       name='Group B: Polar, Uncharged Amino Acids (Hydrophilic)'
   )
   trace3 = go.Bar(
       x=['ASP(D): Aspartic Acid', 'GLU(E): Glutamic Acid'],
       y=[Asp, Glu],
       name='Group C: Polar, Negatively Charged Amino Acids (Hydrophilic)'
   )
   trace4 = go.Bar(
       x=['Lys(K): Lysine','ARG(R): Arganine','HIS(H): Histidine'],
       y=[Lys, Arg,His],
       name='Group C: Polar, Positive Charged Amino Acids (Hydrophilic)'
   )

   data = [trace1, trace2,trace3,trace4]


   ply.plot(data, filename='simple_plot.jpg')


   #

   ListaX=[]
   ListaY = []
   ListaZ=[]

   ca_pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
   for pdb_file in nombrearchivo:
       filename=os.path.basename(pdb_file).split('.')[0]
       chain_dict=dict()
       chain_list=[]

       fp=open(nombrearchivo,'rU')
       for line in fp.read().splitlines():
           if line.startswith("ATOM"):
               
               x1=float(line[32:38])
               y1=float(line[40:46])
               z1=float(line[48:54])
               print(x1)
               print(y1)
               print(z1)
               ListaX.append(x1)
               ListaY.append(y1)
               ListaZ.append(z1)
               
           elif line.startswith("ENDMDL"):
               break
               
       fp.close()



   trace = go.Scatter3d(
       x=list(ListaX), y=list(ListaY), z=list(ListaZ),
       marker=dict(
           size=4,
           color=ListaZ,
           colorscale='Viridis',
       ),
       line=dict(
           color='#1f77b4',
           width=1
       )
   )

   data = [trace]

   ply.plot(data, filename='simple_plot.jpg')
    




root = tk.Tk()


text2 = Label(root, text="Herramienta bioinformatica para analizar archivos PDB")
text2.pack()


text3 = Label(root, text="PDB id:")
text3.pack()

entry_text = tk.StringVar()
entry = tk.Entry(root, width=10, textvariable=entry_text)
entry.pack()



button = tk.Button(root, text="Ejecutar", command=show_text)
button.pack()

label_text = tk.StringVar()
label = tk.Label(root, textvariable=label_text)
label.pack()

root.mainloop()
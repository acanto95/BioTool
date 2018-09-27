#!/usr/bin/env python
import sys,os
import webbrowser, os
import re
import numpy as np
import urllib.request
import plotly.graph_objs as go
import plotly.offline as ply
import tkinter as tk
import PIL
from PIL import Image, ImageDraw, ImageFont
import text_to_image
from tkinter import *
from Bio.PDB import *
from Bio import SeqIO

data = ''


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


   parser = PDBParser()
   structure = parser.get_structure(nombre, nombrearchivo)


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


       datafasta= ''
       datafastaformat = ''

       
       for chain in chain_list:
          datafasta = ('M'+chain + chain_dict[chain])

       datafastaformat = datafasta

       datafastaformat = '<br>'.join(datafastaformat[i:i+80] for i in range(0, len(datafastaformat), 80))





       








   percentlen = len(datafasta)

   Asp = float(datafasta.count('D'))
   Thr = float(datafasta.count('T'))
   Esr = float(datafasta.count('S'))
   Glu = float(datafasta.count('E'))
   Pro = float(datafasta.count('P'))
   Gly = float(datafasta.count('G'))
   Ala = float(datafasta.count('A'))
   Cys = float(datafasta.count('C'))
   Vla=  float(datafasta.count('V'))
   Met = float(datafasta.count('M'))
   Ile = float(datafasta.count('I'))
   Leu = float(datafasta.count('L'))
   Tyr = float(datafasta.count('Y'))
   Phe = float(datafasta.count('F'))
   His = float(datafasta.count('H'))
   Lys = float(datafasta.count('K'))
   Arg = float(datafasta.count('R'))
   Trp = float(datafasta.count('W'))
   Gln = float(datafasta.count('Q'))
   Asn = float(datafasta.count('N'))


   Asp = ((Asp*100)/percentlen)
   Thr = ((Thr*100)/percentlen)
   Esr = ((Esr*100)/percentlen)
   Glu = ((Glu*100)/percentlen)
   Pro = ((Pro*100)/percentlen)
   Gly = ((Gly*100)/percentlen)
   Ala = ((Ala*100)/percentlen)
   Cys = ((Cys*100)/percentlen)
   Vla = ((Vla*100)/percentlen)
   Met = ((Met*100)/percentlen)
   Ile = ((Ile*100)/percentlen)
   Leu = ((Leu*100)/percentlen)
   Tyr = ((Tyr*100)/percentlen)
   Phe = ((Phe*100)/percentlen)
   His = ((His*100)/percentlen)
   Lys = ((Lys*100)/percentlen)
   Arg = ((Arg*100)/percentlen)
   Trp = ((Trp*100)/percentlen)
   Gln = ((Gln*100)/percentlen)
   Asn = ((Asn*100)/percentlen)


   n = 201
   x = np.linspace(0, 2.0*np.pi, n)
   y1 = np.sin(x)
   y2 = np.cos(x)
   y3 = y1 + y2



   ca_pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
   for pdb_file in nombrearchivo:
       filename=os.path.basename(pdb_file).split('.')[0]
       chain_dict=dict()
       chain_list=[]

       fp=open(nombrearchivo,'rU')
       for line in fp.read().splitlines():
           if "MOLECULE" in line:
               
               title = line[10:80]
               
           elif "CHAIN" in line:
               break
               
       fp.close()


   trace1 = go.Bar(
       x=['<a href="http://www.aminoacidsguide.com/Gly.html"> GLY(G): Glycine</a>', '<a href="http://www.aminoacidsguide.com/Ala.html">ALA(A): Alanine</a>','<a href="http://www.aminoacidsguide.com/Val.html"> VAL(V): Valine</a>','<a href="http://www.aminoacidsguide.com/Leu.html">LEU(L): Leucine</a>','<a href="http://www.aminoacidsguide.com/Ile.html">ILE(I): Isoleucine</a>','<a href="http://www.aminoacidsguide.com/Met.html">MET(M): Methionine</a>','<a href="http://www.aminoacidsguide.com/Phe.html">PHE(F): Phenylalanine</a>','<a href="http://www.aminoacidsguide.com/Trp.html">TRP(W): Tryptophan</a>','<a href="http://www.aminoacidsguide.com/Pro.html">PRO(P): Proline</a>'],
       y=[Gly, Ala, Vla, Leu,Ile,Met,Phe,Trp,Pro],
       name='Group A: Nonpolar Amino Acids (Hydrophobic)',
       showlegend=True
   )
   trace2 = go.Bar(
       x=['<a href="http://www.aminoacidsguide.com/Ser.html">SER(S): Serine</a>', '<a href="http://www.aminoacidsguide.com/Thr.html">THR(T): Threonine</a>', '<a href="http://www.aminoacidsguide.com/Cysteine.html">CYS(C): Cysteine</a>','<a href="hhttp://www.aminoacidsguide.com/Tyr.html">TYR(Y): Tyrosine</a>','<a href="http://www.aminoacidsguide.com/Asn.html">ASN(N): Asparagine</a>','<a href="http://www.aminoacidsguide.com/Gln.html">GLN(Q): Glutamine</a>'],
       y=[Esr, Thr, Cys,Tyr, Asn, Gln],
       name='Group B: Polar, Uncharged Amino Acids (Hydrophilic)',
       showlegend=True
   )
   trace3 = go.Bar(
       x=['<a href="http://www.aminoacidsguide.com/Asp.html">ASP(D): Aspartic Acid</a>', '<a href="http://www.aminoacidsguide.com/Glu.html">GLU(E): Glutamic Acid</a>'],
       y=[Asp, Glu],
       name='Group C: Polar, Negatively Charged Amino Acids (Hydrophilic)',
       showlegend=True
   )
   trace4 = go.Bar(
       x=['<a href="http://www.aminoacidsguide.com/Lys.html">Lys(K): Lysine</a>','<a href="http://www.aminoacidsguide.com/Arg.html">ARG(R): Arganine</a>','<a href="http://www.aminoacidsguide.com/His.html">HIS(H): Histidine</a>'],
       y=[Lys, Arg,His],
       name='Group C: Polar, Positive Charged Amino Acids (Hydrophilic)',
       showlegend=True
   )


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


   data = [trace,trace1,trace2,trace3,trace4]
   layout = {
  "plot_bgcolor": 'black',
  "paper_bgcolor": 'black',
  "titlefont": {
      "size": 15,
      "family": "Raleway"
  },
  "font": {
      "color": 'white'
  },
  "margin": {
    "r": 10,
    "t": 20,
    "b": 20,
    "l": 10
  },
  "scene": {"domain": {
      "x": [0, 0.38],
      "y": [0.2, 1]
    },
           "xaxis": {"gridcolor": 'white'},
           "yaxis": {"gridcolor": 'white'},
           "zaxis": {"gridcolor": 'white'}
           },
  "showlegend": True,
  "legend":dict(x=-.1, y=1.2),
  "title": title,
  "xaxis": {
    "anchor": "y",
    "domain": [0.6, 0.93]
  },
  "yaxis": {
    "anchor": "x",
    "domain": [0.2, 1],
    "showgrid": False
  }
   }
   annotations = { "text":"FASTA: " + datafastaformat ,
               "showarrow": False,
               "xref": "paper",
               "yref": "paper",
               "x": -0.06,
               "y": -0.04}

   layout['annotations'] = [annotations]



   fig = go.Figure(data=data, layout=layout)

   ply.plot(fig, filename='simple_plot.jpg')

   encoded_image_path = text_to_image.encode("Hello World!", "image.png")
  



    




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
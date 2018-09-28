#!/usr/bin/env python
import sys,os
import webbrowser, os
import re
import numpy as np
import urllib.request
import plotly.graph_objs as go
import plotly.offline as ply
import tkinter as bioTool
import PIL
from tkinter import *
from Bio.PDB import *
from Bio import SeqIO
from PIL import ImageTk, Image

Estructura={"A":[1.25,0.89,0.78],"R":[0.99,1.02,0.88],"N":[0.87,0.86,1.28],"D":[1.03,0.74,1.41], "C":[1.12,0.85,0.8],"E":[1.45,0.65,1], "Q":[1.24,0.82,0.97], "G":[0.57,0.93,1.64], "H":[1.25,1.04,0.69], "I":[0.94,1.41,0.51], "L":[1.32,1.03,0.59], "K":[1.24,0.81,0.96], "M":[1.43,0.99,0.39], "F":[1.08,1.22,0.58], "P":[0.6,0.71,1.91], "S":[0.82,0.96,1.33], "T":[0.81,1.13,1.03], "W":[1.03,1.15,0.75], "Y":[0.75,1.25,1.05], "V":[0.88,1.48,0.47]}

data = ''


def probAA (Text):

  prob={}
  
  k=3
  n=len(Text)
  
  for i in range (0,n-k+1,3):
  
    pattern=Text[i:i+k]
    
    prob[pattern]=[0,0,0]
    
    prob[pattern][0]=round((Estructura[pattern[0]][0]+Estructura[pattern[1]][0]+Estructura[pattern[2]][0])/3,4)
    
    prob[pattern][1]=round((Estructura[pattern[0]][1]+Estructura[pattern[1]][1]+Estructura[pattern[2]][1])/3,4)
    
    prob[pattern][2]=round((Estructura[pattern[0]][2]+Estructura[pattern[1]][2]+Estructura[pattern[2]][2])/3,4)
    
  return prob 

  


#print (Estructura["A"][0]+Estructura["R"][0])


#list1=[0.98,0.7,1.18]

#criterios de asignación de estructura alfa, beta, girobeta y azar

def estrAA (Text):

  prob=probAA(Text)
  
  estr={}
  
  for pattern in prob:
    if prob[pattern][0]>1.1 and prob[pattern][1]<1.2 and prob[pattern][2]<1.3:
      estr[pattern]="alfa"
    else:
      if prob[pattern][0]<1.1 and prob[pattern][1]>1 and prob[pattern][2]<1.3:
        estr[pattern]="beta"
      else:
        if prob[pattern][0]<1.25 and prob[pattern][1]<1 and prob[pattern][2]>1.15:
          estr[pattern]="giro beta"
        else:
          estr[pattern]="azar"
  return estr
  
#print(estrAA(list1))

#regresa los patrones ordenados por estructura

def countEstr(Text):
  alfa=[]
  
  beta=[]
  
  bgiro=[]
  
  azar=[]
  
  estr=estrAA(Text)
  
  for pattern in estr:
  
    if estr[pattern]=="alfa":
    
      alfa.append(pattern)
      
    else:
    
      if estr[pattern]=="beta":
      
        beta.append(pattern)
        
      else :
      
        if estr[pattern]=="giro beta":
        
          bgiro.append(pattern)
          
        else:
        
          azar.append(pattern)

          listaresult = [alfa,len(alfa),beta,len(beta),bgiro,len(bgiro),azar,len(azar)]
          
  return listaresult
  
#regresa las listas y su cardinalidad



def build_guipro():    #Funcion principal que es llamada en la interfaz grafica


   nombre= entry_text.get() # Recibe datos de la GUI
   nombrearchivo = nombre + ".pdb"
   connection=urllib.request.urlopen("https://files.rcsb.org/view/"+nombre+".pdb") #Descarga archivo PDB
   output=connection.read().decode('utf-8')


   # Escribit datos al archivo
   filename = nombrearchivo
   file2 = open(filename, 'w+')
   file2.write(str(output))
   file2.close()



   aa3to1={   #Diccionario para Aminoacidos
      'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
      'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
      'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
      'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
      'MSE':'M',
   }
                #Guia para el lector de archivos para que sepa donde buscar 
   ca_pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
   for pdb_file in nombrearchivo: 
       filename=os.path.basename(pdb_file).split('.')[0]
       chain_dict=dict()
       chain_list=[]

       #Concatenacion de valores PDB a FASTA
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
       datafastaformat = '' ##Variables para darle formato diferente al FASTA

       #Se agrega el Aminoacido "M" a la secuencia y se agrega a un string
       for chain in chain_list:
          datafasta = ('M'+chain + chain_dict[chain])

       datafastaformat = datafasta

       #Transformacion del string FASTA para que la GUI lo pueda leer 
       datafastaformat = '<br>'.join(datafastaformat[i:i+80] for i in range(0, len(datafastaformat), 80))





       






   #Calculo de porcentajes

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


   #Lectura de la molecula para la GUI
   ca_pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
   for pdb_file in nombrearchivo:
       filename=os.path.basename(pdb_file).split('.')[0]
       chain_dict=dict()
       chain_list=[]

       fp=open(nombrearchivo,'rU')
       for line in fp.read().splitlines():
           if "MOLECULE" in line: #quitar nombre molecule
               
               title = line[21:80]
               
           elif "CHAIN" in line:
               break
               
       fp.close()

       title = title + "<br>" + "Numero total de aminoácidos: "+ str(percentlen)

       #Graficas
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
   #Grafica en 3D
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
      "x": [0.52, 0.97],
      "y": [0.3, 1]
    },
           "xaxis": {"gridcolor": 'white'},
           "yaxis": {"gridcolor": 'white'},
           "zaxis": {"gridcolor": 'white'}
           },
  "showlegend": True,
  "legend":dict(x=0, y=1.2),
  "title": title,
  "xaxis": {
    "anchor": "y",
    "domain": [0.01, 0.45]
  },
  "yaxis": {
    "anchor": "x",
    "domain": [0.26, .95],
    "showgrid": False
  }
   }
   annotations = { "text":"FASTA: " + datafastaformat +'<br>'+'<a href="http://www.bachem.com/fileadmin/user_upload/pdf/Flyers/Periodic_Chart_Amino_Acids.pdf"><b>Tablas de Aminoácidos</b></a>' ,
               "showarrow": False,
               "xref": "paper",
               "yref": "paper",
               "x": 1,
               "y": 0}

   layout['annotations'] = [annotations]



   fig = go.Figure(data=data, layout=layout)

   ply.plot(fig, filename='simple_plot.jpg')



   lista=countEstr(datafasta)
   mensaje1=" ".join(lista[0])
   mensaje1="Alfa: "+mensaje1
   mensaje2=" ".join(lista[2])
   mensaje2="Beta: "+mensaje2
   mensaje3=" ".join(lista[4])
   mensaje3="Beta giro: "+mensaje3
   mensaje4=" ".join(lista[6])
   mensaje4="Azar : "+mensaje4
      
   labels = [mensaje1,mensaje2,mensaje3,mensaje4]
   values = [lista[1],lista[3],lista[5],lista[7]]
   trace = go.Pie(labels=labels, values=values)
   ply.plot([trace], filename='basic_pie_chart')



  



    



#-----------------------------------Interfaz Grafica
root = bioTool.Tk()


img = ImageTk.PhotoImage(Image.open("3.png"))
panel = bioTool.Label(root, image = img)
panel.pack(side = "top", fill = "both", expand = "yes")

text2 = Label(root, text="Herramienta bioinformática para analizar proteínas")
text2.pack()
text3 = Label(root, text="PDB ID:")
text3.pack()
entry_text = bioTool.StringVar()
entry = bioTool.Entry(root, width=10, textvariable=entry_text)
entry.pack()





button = bioTool.Button(root, text="Desplegar Información", command=build_guipro, height = 3, width = 35)
button.pack()

label_text = bioTool.StringVar()
label = bioTool.Label(root, textvariable=label_text)
label.pack()

root.title("BioTool")


root.mainloop()
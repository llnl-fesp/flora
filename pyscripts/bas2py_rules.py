from flora import *
subrules = [
   ['\(','['],
   ['\)',']'],
   #[';','\n'],
   ['^!','#!'],
   ['^ *',''],
   ['^\t*',''],
   [r'\bgenerate\b','glr.glrgen()'],
   [r'\btrue\b','True'],
   [r'\bfalse\b','False'],
   ]
warnrules = []
def raw_string(s):
    s = s.encode('unicode-escape').decode()
    return s

for p in [('glr',flora.glr),('utl',flora.utl)]:
   for v in p[1].varlist():
       subrules.append([r'\b'+raw_string(v)+r'\b',p[0]+'.'+v])

       if "Dimension:" in p[1].listvar(v):
          d = p[1].listvar(v).split("Dimension:")[1].split("\n")
          if "0:" in d[0]:
             warnrules.append([r'\b'+raw_string(v)+r'\b','base 0, '+d[0]])

       


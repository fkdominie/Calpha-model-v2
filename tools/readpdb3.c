/****************************************************************************/
/* readpdb.c                                                                */
/* flag = 1, Ca atoms                                                       */
/* flag = 2, N,Ca,C atoms                                                   */
/* flag = 3, all atoms (N,Ca,C,O + non-H sidechain atoms)                   */
/****************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# define NMAX 200 
# define S1 " N  "," H  "," CA "," C  "," O  "
# define S2 " N  "," H  "," CA "," C  "," O  "," CB "

const char atom[20][15][5]={
{S1},                                                                /* gly */
{S2},                                                                /* ala */
{S2," CG1"," CG2"},                                                  /* val */
{S2," CG "," CD1"," CD2"},                                           /* leu */
{S2," CG1"," CG2"," CD1"},                                           /* ile */
{S2," OG "},                                                         /* ser */
{S2," OG1"," CG2"},                                                  /* thr */
{S2," SG "},                                                         /* cys */
{S2," CG "," SD "," CE "},                                           /* met */
{S2," CG "," CD "},                                                  /* pro */
{S2," CG "," OD1"," OD2"},                                           /* asp */
{S2," CG "," OD1"," ND2"},                                           /* asn */
{S2," CG "," CD "," OE1"," OE2"},                                    /* glu */
{S2," CG "," CD "," OE1"," NE2"},                                    /* gln */
{S2," CG "," CD "," CE "," NZ "},                                    /* lys */
{S2," CG "," CD "," NE "," CZ "," NH1"," NH2"},                      /* arg */
{S2," CG "," ND1"," CD2"," CE1"," NE2"},                             /* his */
{S2," CG "," CD1"," CD2"," CE1"," CE2"," CZ "},                      /* phe */
{S2," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH "},               /* tyr */
{S2," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"}  /* trp */
};

const char amino[20][4]={"GLY","ALA","VAL","LEU","ILE","SER","THR",
			 "CYS","MET","PRO","ASP","ASN","GLU","GLN",
			 "LYS","ARG","HIS","PHE","TYR","TRP"};

const char aalett[20]={'G','A','V','L','I','S','T',
		       'C','M','P','D','N','E','Q',
		       'K','R','H','F','Y','W'};

//const int natom[20]={4,5,7,8,8,6,7,6,8,7,8,8,9,9,9,11,10,11,12,14};
const int natom[20]={5,6,8,9,9,7,8,7,9,8,9,9,10,10,10,12,11,12,13,15};

/****************************************************************************/
int readpdb(int iflag,char *fn,double x[],double y[],double z[],int seq[],int aan[]);
void dumppdb(char *fn,double x[],double y[],double z[],int seq[],int n);
/****************************************************************************/
int fgetline(char *line,int max,FILE *fp);
void substr(char *sub,char *str,int pos, int offset);
int get_iaa(char c);
/****************************************************************************/
int readpdb(int iflag,char *fn,double x[],double y[],double z[],int seq[],int aan[]){
  char line[100],str1[20],str2[20],str3[20],str4[20];
  int i0,aa,atm,iaa=0,iatm; 
  FILE *fp; 

  if (iflag==0) return 0;

  fp=fopen(fn,"r");

  /* read through until first ATOM block */

  do {
    fgetline(line,100,fp);
    substr(str1,line,0,6);
    substr(str2,line,12,4); // atom name 
    substr(str3,line,17,3); // amino acid name 
    substr(str4,line,22,4); // amino acid number 
  } while ( strcmp(str1,"ATOM  ") != 0);


  aa = -1;
  atm = 0;
  i0 = 0;
  
  while ( strcmp(str1,"ATOM  ") == 0 && feof(fp) != EOF ) {
      
    /* found next amino acid? */

    if (strcmp(str2," N  ") == 0) {    
      if (aa > -1 && atm != natom[iaa] && iflag == 3) {
	printf("# !! missing atom(s) at position %i %s \n",aa,amino[iaa]);
	atm = natom[iaa];
      }
      
      aa++;             // aa counter
      i0 += atm;        // atom counter (total)
      atm = 0;          // reset atom counter (current aa)

      /* which amino acid type was found? */
      
      iaa=0;  while ( strcmp(amino[iaa],str3)!=0 && iaa<20) ++iaa;

      if (iaa >= 20) {
	printf("readpdb -- unknown amino acid\n");
	fclose(fp);
	exit(-1);
      }

      seq[aa] = aalett[iaa];
      aan[aa] = atoi(str4);
    }

    /* which atom type was found? */
    
    iatm=0;  while ( strcmp(str2,atom[iaa][iatm])!=0 && iatm<natom[iaa]) iatm++; 
    
    if (iatm >= natom[iaa]) {
      //      printf("Unknown atom, amino acid %s %d %c \n",str2,aa,seq[aa]);

      if( fgetline(line,100,fp) == 0) {fclose(fp); return aa+1;}
      substr(str1,line,0,6);
      substr(str2,line,12,4);      
      substr(str3,line,17,3);
      substr(str4,line,22,4); 
      
      continue;
    }
    
    /* only CA atoms (iflag = 0) */
    
    if ( (iatm == 1 && iflag == 1) ) {
      substr(str4,line,31,8); x[aa] = atof(str4);
      substr(str4,line,39,8); y[aa] = atof(str4);
      substr(str4,line,47,8); z[aa] = atof(str4);
      //      printf("%s %f %f %f %i\n",amino[iaa],x[aa],y[aa],z[aa],aan[aa]);
    }
    
    /* N/CA/C atoms (iflag=2) or All atoms (iflag=3) */

    if ( (iatm >= 0 && iatm <= 2 && iflag == 2) || (iflag == 3) ) {
      substr(str4,line,31,8); x[i0 + iatm] = atof(str4);
      substr(str4,line,39,8); y[i0 + iatm] = atof(str4);
      substr(str4,line,47,8); z[i0 + iatm] = atof(str4);
      //      printf("%s %f %f %f %i\n",amino[iaa],x[i0+iatm],y[i0+iatm],z[i0+iatm],aan[aa]);
      atm++;
    }

    /* read next line */
    
    if (fgetline(line,100,fp) == 0) {fclose(fp); return aa + 1;}
    substr(str1,line,0,6);
    substr(str2,line,12,4);      
    substr(str3,line,17,3);
    substr(str4,line,22,4); 
  }

  fclose(fp);

  return aa + 1;
}
/****************************************************************************/
void substr(char *sub,char *str,int pos, int offset){
  strncpy(sub,str+pos,offset);
  *(sub+offset)='\0';
}
/****************************************************************************/
int fgetline(char *line,int max,FILE *fp){
  if (fgets(line,max,fp) == NULL){
    strcpy(line,"");
    return 0;
  }
  else
    return strlen(line);
}
/****************************************************************************/
int get_iaa(char c) {
  int i;
  for (i=0;i<20;i++) if (aalett[i] == c) return i;
  return -1;
}
/****************************************************************************/
void dumppdb(char *fn,double x[],double y[],double z[],int seq[],int n) {
  int i,j,k=0,iaa;
  FILE *fp;

  fp=fopen(fn,"w");

  for (i=0;i<n;i++) {
    iaa=get_iaa(seq[i]);
    for (j=0;j<natom[iaa];j++) {
      fprintf(fp,"ATOM  %4u  %s %s  %4u    %8.3f%8.3f%8.3f  1.00\n",
	      k+1,atom[iaa][j],amino[iaa],i+1,x[k],y[k],z[k]); 
      k++;
    }
  }

  fclose(fp);
}
/****************************************************************************/

#undef S1 
#undef S2












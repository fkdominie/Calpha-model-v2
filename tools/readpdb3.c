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

# define S1 " N  "," CA "," C  "," O  "
# define S2 " N  "," CA "," C  "," O  "," CB "

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

const int natom[20]={4,5,7,8,8,6,7,6,8,7,8,8,9,9,9,11,10,11,12,14};

/****************************************************************************/
int readpdb(int iflag,char *fn,double x[],double y[],double z[],int seq[]);
void dumppdb(char *fn,double x[],double y[],double z[],int seq[],int n);
/****************************************************************************/
int fgetline(char *line,int max,FILE *fp);
void substr(char *sub,char *str,int pos,int len);
int get_iaa(char c);
int get_iaa3(char *str);
int get_iatm(char *str,int iaa);
/****************************************************************************/
int readpdb(int iflag,char *fn,double x[],double y[],double z[],int seq[]){
  char line[100],str1[20],str2[20],str3[20],str4[20];
  int i0,aa,atm,iaa=0,iatm; 
  FILE *fp; 

  /* iflag == 1 : read CA atoms
     iflag == 3 : read all non-H atoms */

  if (iflag == 0)
    return 0;

  fp = fopen(fn,"r");

  /* Find the ATOM block */

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
      
    if (strcmp(str2," N  ") == 0) {    

      if (aa > -1 && atm != natom[iaa] && iflag == 3) {
	printf("# !! missing atom(s) at position %i %s \n",aa,amino[iaa]);
	atm = natom[iaa];
      }
      
      aa++;             // aa counter
      i0 += atm;        // atom counter
      atm = 0;          // reset

      if ( (iaa = get_iaa3(str3)) == 20) {
	printf("readpdb -- unknown amino acid at position %i\n",aa);
	fclose(fp);
	exit(-1);
      }
      
      seq[aa] = aalett[iaa];
    }
    
    if ( (iatm = get_iatm(str2,iaa)) == natom[iaa] )  {
      
      if ( fgetline(line,100,fp) == 0 ) {
	fclose(fp);
	return aa + 1;
      }
      substr(str1,line,0,6);
      substr(str2,line,12,4);      
      substr(str3,line,17,3);
      substr(str4,line,22,4); 
      continue;
    }
    
    if ( iatm == 1 && iflag == 1 ) {  /* CA atoms */
      substr(str4,line,31,8); x[aa] = atof(str4);
      substr(str4,line,39,8); y[aa] = atof(str4);
      substr(str4,line,47,8); z[aa] = atof(str4);
      //      printf("%s %f %f %f %i\n",amino[iaa],x[aa],y[aa],z[aa],aan[aa]);
    }
    
    if ( iflag == 3 ) { /* non-H atoms */
      substr(str4,line,31,8); x[i0 + iatm] = atof(str4);
      substr(str4,line,39,8); y[i0 + iatm] = atof(str4);
      substr(str4,line,47,8); z[i0 + iatm] = atof(str4);
      //      printf("%i %i %i %s %f %f %f\n",iaa,i0+iatm,atm,amino[iaa],x[i0+iatm],y[i0+iatm],z[i0+iatm]);
      atm++;
    }
    
    if (fgetline(line,100,fp) == 0) {     /* read next line */
      fclose(fp);
      return aa + 1;
    }
    substr(str1,line,0,6);
    substr(str2,line,12,4);      
    substr(str3,line,17,3);
    substr(str4,line,22,4); 
  }

  fclose(fp);

  return aa + 1;
}
/****************************************************************************/
void dumppdb(char *fn,double x[],double y[],double z[],int seq[],int n) {
  int i,j,k = 0,iaa;
  FILE *fp;

  fp = fopen(fn,"w");

  for (i = 0; i < n; i++) {
    iaa = get_iaa(seq[i]);
    for (j = 0; j < natom[iaa]; j++) {
      fprintf(fp,"ATOM  %4u  %s %s  %4u    %8.3f%8.3f%8.3f  1.00\n",
	      k+1,atom[iaa][j],amino[iaa],i+1,x[k],y[k],z[k]); 
      k++;
    }
  }

  fclose(fp);
}
/****************************************************************************/
void substr(char *sub,char *str,int pos,int len){
  int i;
  for (i = 0; i < len; i++) {
    if ( (sub[i] = str[pos + i]) == '\0')
      return;
  }
  sub[len] = '\0';
}
/****************************************************************************/
int fgetline(char *line,int max,FILE *fp){
  if (fgets(line,max,fp) == NULL) {
    strcpy(line,"");
    return 0;
  } 
  return strlen(line);
}
/****************************************************************************/
int get_iaa(char c) {
  int i;
  for (i = 0; aalett[i] != c && i < 20; ++i);
  return i;
}
/****************************************************************************/
int get_iaa3(char *str) {
  int i;
  for (i = 0; strcmp(amino[i],str) != 0 && i < 20; ++i);
  return i;
}
/****************************************************************************/
 int get_iatm(char *str,int iaa) {
  int i;
  for (i = 0; strcmp(str,atom[iaa][i]) != 0 && i < natom[iaa]; ++i);
  return i;
}
/****************************************************************************/

#undef S1 
#undef S2












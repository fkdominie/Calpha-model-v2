# include <stdio.h>
# include <stdio.h>

int i,j;
int cc[500][500];
FILE *fp;

int main(int argc, char *argv[]) {

  fp = fopen("smog_2jp1_9-51","r");
  while (2 == fscanf(fp,"%i %i",&i,&j))
    cc[i][j] = 1;
  fclose(fp);

  for (i = 0; i < 93; i++) {
    for (j = i+1; j < 93; j++) {
      if (cc[i][j] == 0 && cc[i+93][j+93] == 1)
	printf("%i %i\n",i,j);
    }
  }

  for (i = 93; i < 186; i++) {
    for (j = i+1; j < 186; j++) {
      if (cc[i][j] == 0 && cc[i-93][j-93] == 1)
	printf("%i %i\n",i,j);
    }
  }

  
  /*  for (i = 0; i < 93; i++) {
    for (j = i+1; j < 93; j++) {
      if (cc[i][j] || cc[i+93][j+93])
	printf("%i %i\n",i,j);
    }
  }

  for (i = 93; i < 186; i++) {
    for (j = i+1; j < 186; j++) {
      if (cc[i][j] || cc[i-93][j-93])
	printf("%i %i\n",i,j);
    }
  }

  for (i = 0; i < 93; i++) {  
    for (j = 93; j < 186; j++) {
      if (cc[i][j])
	printf("%i %i\n",i,j);
    }
    } */
  
  return 0;
}

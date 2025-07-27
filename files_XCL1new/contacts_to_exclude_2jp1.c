# include <stdio.h>
# include <stdio.h>

int i,j;
int cc[500][500];
FILE *fp;

int main(int argc, char *argv[]) {

  fp = fopen("smog_2jp1_r8-52","r");
  while (2 == fscanf(fp,"%i %i",&i,&j))
    cc[i][j] = 1;

  for (i = 0; i < 93; i++) {
    for (j = i+1; j < 93; j++) {
      if (cc[i][j] == 1 && cc[i+93][j+93] == 0)
	printf("(%i,%i); ",i,j);
    }
  }

  for (i = 93; i < 186; i++) {
    for (j = i+1; j < 186; j++) {
      if (cc[i][j] == 1 && cc[i-93][j-93] == 0)
	printf("(%i,%i); ",i,j);
    }
  }
  printf("\n");

  return 0;
  
  for (i = 0; i < 93; i++) {
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
  }
  
  fclose(fp);
  return 0;
}

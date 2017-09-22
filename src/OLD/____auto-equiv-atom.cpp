#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ATOM	(1024)

void Quit_With_Error_Msg(char szMsg[]);

int nAtom;

int main(void)
{
	FILE *fIn;
	int i, j, Representative[MAX_ATOM], Output[MAX_ATOM], nAtom, Index, ReadItem;

	fIn = fopen("equiv.txt", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file equiv.txt.\nQuit\n");
	}

	memset(Output, 0, sizeof(int)*MAX_ATOM);

	nAtom = 0;
	while(1)	{
		ReadItem = fscanf(fIn, "%d %d", &Index, &(Representative[nAtom]));
		if(ReadItem == 2)	{
			if(Representative[nAtom] > 0)	{	// a valid new entry for atom equivalency
				Representative[nAtom]--;
				Output[Representative[nAtom]] = 1;
			}
			else	{
				Representative[nAtom] = -1;
			}
			nAtom++;
		}
		else	{
			break;
		}
	}
	fclose(fIn);

	for(i=0; i<nAtom; i++)	{
		if(Output[i] == 1)	{
			printf("equivalent  %3d", i+1);
			for(j=0; j<nAtom; j++)	{
				if(j==i)	{
					continue;
				}
				if(Representative[j] == i)	{
					printf(" %3d", j+1);
				}
			}
			printf("\n");
		}
	}


	return 0;
}


void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("../error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in auto-equiv-atom.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}



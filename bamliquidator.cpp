#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include "sam.h"

/* The MIT License (MIT) 

   Copyright (c) 2013 Xin Zhong and Charles Lin

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE. 
 */

struct genericItem
{
  /* to include beditem and readitem 
  when only coordinate is concerned
  and when bed/sam format cannot be sured
  */
  struct genericItem *next;
  unsigned int start;
  unsigned int stop;
};


struct readitem
{
  struct readitem *next;
  unsigned int start;
  /* read stop is start + strlen(seq)
  this *stop* will only be used for computing density
  will not be reported to js for bed plotting
  the actual stop need to be determined by cigar
  */
  unsigned int stop;
  char *id;
  uint32_t flag; // flag from bam
  char strand;
  char *flag_str; // flag str from sam
  char *seq; // read sequence
  uint32_t *cigar; // cigar from bam
  uint32_t n_cigar;
  char *cigar_str; // cigar str from sam
  char *mismatch;
};


int intMin(int a, int b)
{
  if(a < b) return a;
  return b;
}

int intMax(int a, int b)
{
  if(a > b) return a;
  return b;
}



struct userData
{
  struct readitem *sl;
  char strand;
  unsigned int extendlen;
};




static int bam_fetch_func(const bam1_t *b,void *data)
{
  if (b->core.tid < 0) return 0;

  struct userData *udata=(struct userData *)data;

  const bam1_core_t *c = &b->core;

  char strand= (c->flag&BAM_FREVERSE)?'-':'+';
  if(udata->strand=='+')
  {
    if(strand!='+') return 0;
  }
  else if(udata->strand=='-')
  {
    if(strand!='-') return 0;
  }

  struct readitem *r=(struct readitem*) malloc(sizeof(struct readitem));
  r->strand=strand;

  r->cigar_str=NULL; // sam
  r->n_cigar=c->n_cigar;
  r->cigar=(uint32_t*) malloc(sizeof(uint32_t)*c->n_cigar);
  uint32_t *cigar = bam1_cigar(b);

  // get read length
  int i, readlen;
  if (b->core.tid < 0) return 0;
  for (i = readlen = 0; i < c->n_cigar; ++i)
  {
    int op = cigar[i]&0xf;
    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
      readlen += cigar[i]>>4;
  }



  for (i=0; i < c->n_cigar; ++i)
  {
    (r->cigar)[i]=cigar[i]; // copy cigar
  }

  uint8_t *seq=bam1_seq(b);
  r->seq=(char*) malloc(sizeof(char)*(c->l_qseq+1));
  for(i=0; i<c->l_qseq; i++)
  {
    (r->seq)[i]=bam_nt16_rev_table[bam1_seqi(seq,i)];
  }
  (r->seq)[i]='\0';

  r->id=strdup(bam1_qname(b));
  r->start=c->pos;
  r->stop=c->pos+readlen;

  //printf("%d\t%d\t%c\t", r->start, r->stop, strand);

  // extend
  if(udata->extendlen>0)
  {
    if(strand=='+')
    {
      r->stop+=udata->extendlen;
    }
    else
    {
      r->start=intMax(0,r->start-udata->extendlen);
    }
  }

  //printf("%d\t%d\n", r->start, r->stop);

  r->flag=c->flag;
  r->flag_str=NULL;
  uint8_t *m=bam_aux_get(b,"MD");
  if(m)
  {
    r->mismatch=strdup((char *)m);
  }
  else
  {
    r->mismatch=NULL;
  }
  r->next=udata->sl;
  udata->sl=r;
  return 0;
}



struct readitem *bamQuery_region(samfile_t *fp, bam_index_t *idx, char *coord, char strand, unsigned int extendlen)
{
  // will not fill chromidx
  int ref,beg,end;
  int rc = bam_parse_region(fp->header,coord,&ref,&beg,&end);
  if (rc != 0)
  {
    printf("bam_parse_region failed with return code %d, exiting\n", rc);
    exit(rc);
  }
  if(ref<0)
  {
    return NULL;
  }
  struct userData *d=(struct userData*) malloc(sizeof(struct userData));
  d->sl=NULL;
  d->strand=strand;
  d->extendlen=extendlen;
  bam_fetch(fp->x.bam,idx,ref,beg,end,d,bam_fetch_func);
  return d->sl;
}


char *getDepositePath4url(char *url)
{
  char *urlcopy=strdup(url);
  char delim[]="/";
  char *tok=strtok(urlcopy,delim);
  char *deposit_dir=strdup("./");
  struct stat sbuf;
  while(1)
  {
    assert(asprintf(&deposit_dir,"%s/%s",deposit_dir,tok)>0);
    if(stat(deposit_dir,&sbuf)==-1)
    {
      // create new dir
      if(mkdir(deposit_dir,S_IRWXU)!=0)
      {
        fprintf(stderr,"failed to create caching directory\n");
        return NULL;
      }
    }
    tok=strtok(NULL,delim);
    if(tok==NULL)
    {
      break;
    }
  }
  return deposit_dir;
}

/**************** INIT
*/

int main(int argc, char *argv[])
{
  if(argc!=8)
  {
    printf("[ bamliquidator ] output to stdout\n1. bam file (.bai file has to be at same location)\n2. chromosome\n3. start\n4. stop\n5. strand +/-, use dot (.) for both strands\n6. number of summary points\n7. extension length\n\n");
    return 1;
  }

  char *tail=NULL;
  unsigned int start=strtol(argv[3],&tail,10);
  if(tail[0]!='\0')
  {
    fprintf(stderr, "wrong start (%s)\n", argv[3]);
    return 1;
  }
  unsigned int stop=strtol(argv[4],&tail,10);
  if(tail[0]!='\0' || stop<=start)
  {
    fprintf(stderr, "wrong stop (%s)\n", argv[4]);
    return 1;
  }
  char strand=argv[5][0];
  if(strand!='+' && strand!='-' && strand!='.')
  {
    fputs("wrong strand, must be +/-/.\n",stderr);
    return 1;
  }
  unsigned int spnum=strtol(argv[6],&tail,10);
  if(tail[0]!='\0' || spnum<=0)
  {
    fprintf(stderr, "wrong spnum (%s)\n", argv[6]);
    return 1;
  }

  unsigned int extendlen=(unsigned short)strtol(argv[7],&tail,10);
  if(tail[0]!='\0')
  {
    fprintf(stderr, "wrong extension length (%s)\n", argv[7]);
    return 1;
  }

  double *data=(double*) malloc(sizeof(double)*spnum);
  if(data==NULL)
  {
    fputs("out of mem\n", stderr);
    return 1;
  }

  int i;
  for(i=0; i<spnum; i++) data[i]=0;

  samfile_t *fp=NULL;
  bam_index_t *bamidx=NULL;

  char *bamfile=argv[1];

  if((fp=samopen(bamfile,"rb",0))==0)
  {
    printf("samopen() error with %s\n",bamfile);
    return 1;
  }

  if(strncasecmp(bamfile,"http",4)==0)
  {
    char *deposit_dir=getDepositePath4url(bamfile);
    if(chdir(deposit_dir)!=0)
    {
      printf("failed to change to bam caching dir: %s\n", deposit_dir);
      return 1;
    }
  }
  bamidx=bam_index_load(bamfile);

  /* fetch bed items for a region and compute density
  only deal with coord, so use generic item
  */
  int startArr[spnum], stopArr[spnum];
  float pieceLength = (float)(stop-start) / spnum;
  for(i=0; i<spnum; i++)
  {
    startArr[i] = (int)(start + pieceLength*i);
    stopArr[i] = (int)(start + pieceLength*(i+1));
  }

  struct genericItem *itemsl, *item;
  char *tmp;
  if(asprintf(&tmp,"%s:%d-%d",argv[2],start,stop)<0)
  {
    printf("asprintf() out of mem\n");
    return 1;
  }
  itemsl=(struct genericItem *)bamQuery_region(fp,bamidx,tmp,strand,extendlen);
  free(tmp);

  for(item=itemsl; item!=NULL; item=item->next)
  {
    // collapse this bed item onto the density counter
    for(i=0; i<spnum; i++)
    {
      if(item->start > stopArr[i]) continue;
      if(item->stop < startArr[i]) break;
      int start=intMax(item->start,startArr[i]);
      int stop=intMin(item->stop,stopArr[i]);
      if(start<stop)
      {
        // as Charles suggested, add the fraction of the read (overlapping with the bin)
        // instead of just counting the read
        data[i] += stop-start;
      }
    }
  }

  for(i=0; i<spnum; i++)
  {
    printf("%d\n", (int)(data[i]));
  }

  return 0;
}

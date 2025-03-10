#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"

/* here are our X variables */
Display *dis;

int screen;
unsigned char *x_buffer;
Window win;
GC gc;
XImage *x_image;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();

struct tangle { double hyp[2];  double opp[2]; double adj[2]; };


void disp (unsigned char *,int,int);

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}

void plot_tri( char *image, struct tangle *t,int r, int g, int b, int thick)
{

	lineu( image, t->hyp, t->opp,r, g,b,thick);
	lineu( image, t->opp, t->adj,r, g,b,thick);
	lineu( image, t->adj, t->hyp,r, g,b,thick);
	
}

void perp_tri( char *image, struct tangle *t,struct tangle *o1, struct tangle *o2, int r, int g, int b, int thick)
{
	/* Scary maths
	perp slope =-1/m
	Hypotenuse y=mx+c     m=dy/dx
	c=y-mx;
	Perp yp=(-xp/m)+d
       	d=yp+(xp/m);
       */

	double dx,dy,c,d,m,xi,yi,bg[2];

	dx=t->opp[0]-t->hyp[0]; if (dx==0){ dx=0.00001;} // we dont want divide by zero.
	dy=t->opp[1]-t->hyp[1]; if (dy==0){ dy=0.00001;} // we dont want 1/m to be infinity either
	m=dy/dx; c=t->hyp[1]-(m*t->hyp[0]);

	/* The two lines intersect at xn yn perp is at the meeting pount of o and a */
	
	d=t->adj[1]+(t->adj[0]/m);

	/* xi intersect m*xi+c=(-xi/m)+d
	m*xi+(xi/m)=d-c
	xi=(d-c)/(m+(1/m)) */

	xi=(d-c)/(m+(1/m));
	yi=(m*xi)+c;

	bg[0]=xi;bg[1]=yi;
	
	lineu( image, t->adj,bg,r, g,b,thick);

	o1->hyp[0]=t->adj[0];	o1->hyp[1]=t->adj[1];	
	o2->hyp[0]=t->adj[0];	o2->hyp[1]=t->adj[1];	

	o1->opp[0]=t->hyp[0];	o1->opp[1]=t->hyp[1];	
	o2->opp[0]=t->opp[0];	o2->opp[1]=t->opp[1];	

	o1->adj[0]=xi;	o1->adj[1]=yi;	
	o2->adj[0]=xi;	o2->adj[1]=yi;	
}


int qlott(unsigned char *image, int xp, int yp, int r, int g, int b, int xs, int ys)
{
	int x,y;
	if (xs<1){ xs=1;}
	if (ys<1){ ys=1;}
	for (x=0;x<xs;x++)
	{
		for (y=0;y<ys;y++)
		{
			plott(image,xp+x,yp+y,r,g,b);
		}
	}
}	


int main(int argc,char *argv[])
{
	unsigned char *image2,*image3,*image4,*image5,*image6;
	int i,loop;
	struct tangle *tri[500000];
	struct tangle *sri[500000];

	for (i=0;i<500000;i++){
        	tri[i]=(struct tangle *)malloc(sizeof (struct tangle)); // disp buffer
        	sri[i]=(struct tangle *)malloc(sizeof (struct tangle)); // disp buffer
	}



	init_x();

        image2=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        image3=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        image4=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        image5=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
        image6=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer


	int dir,pc;
	char name[100];


	for (dir=0;dir<X_SIZE;dir+=1)
	{	
	clear(image5,0);

	float t,pp;

	t=8*(float)(dir)/X_SIZE;

	pp=(M_PI*sin(M_PI*t/4))/4;


	tri[0]->hyp[0]=(X_SIZE/2)-(Y_SIZE*(sin(2*M_PI*t)/2)); 	tri[0]->hyp[1]=(Y_SIZE/2)-(Y_SIZE*(cos(2*M_PI*t)/2)) ;
	tri[0]->opp[0]=(X_SIZE/2)-(Y_SIZE*(sin((2*M_PI*t)+(2*M_PI/4))/2)); tri[0]->opp[1]=(Y_SIZE/2)-(Y_SIZE*(cos(2*M_PI*t+(2*M_PI/4)))/2) ;
	tri[0]->adj[0]=(X_SIZE/2)-(Y_SIZE*(sin(pp+(2*M_PI*t)+(4*M_PI/4))/2)); tri[0]->adj[1]=(Y_SIZE/2)-(Y_SIZE*(cos(2*M_PI*t+(4*M_PI/4)+pp))/2) ;

	sri[0]->hyp[0]=(X_SIZE/2)-(Y_SIZE*(sin((2*M_PI*t)+(6*M_PI/4))/2)); sri[0]->hyp[1]=(Y_SIZE/2)-(Y_SIZE*(cos(2*M_PI*t+(6*M_PI/4)))/2) ;
	sri[0]->opp[0]=(X_SIZE/2)-(Y_SIZE*(sin((2*M_PI*t)+(4*M_PI/4))/2)); sri[0]->opp[1]=(Y_SIZE/2)-(Y_SIZE*(cos(2*M_PI*t+(4*M_PI/4)))/2) ;
	sri[0]->adj[0]=(X_SIZE/2)-(Y_SIZE*(sin(-pp+(2*M_PI*t)+(0*M_PI/4))/2)); sri[0]->adj[1]=(Y_SIZE/2)-(Y_SIZE*(cos(2*M_PI*t+(0*M_PI/4)-pp))/2) ; 

	int gen;
	gen=0;	
	plot_tri( image5, tri[0],255, 255, 255, 2);
	plot_tri( image5, sri[0],255, 255, 255, 2);
	//perp_tri( image5, tri[gen],tri[gen+1],tri[gen+2],255, 255, 255, 20);
	int num;
	int cc;
	cc=0;
	for (num=0;num<7-(6*sin(2*M_PI*t));num++)
	{
		int p,pn;
		p=pow(2,num);

		for (gen=0;gen<p;gen++)
		{
			//printf("g1 %d g2  %d g3 %d\n",cc,(cc*2)+1,((cc*2)+2));
			perp_tri( image5, tri[cc],tri[(cc*2)+1],tri[(cc*2)+2],128+(127*(sin(t*sqrt(cc)))), 128+(127*(sin(t*(sqrt(2*cc))))),128+(127*(sin(t*sqrt(3*cc)))), 2);
			perp_tri( image5, sri[cc],sri[(cc*2)+1],sri[(cc*2)+2],128+(127*(sin(t*sqrt(cc*2)))), 128+(127*(sin(t*(sqrt(5*cc))))),128+(127*(sin(t*sqrt(7*cc)))), 2);
			cc++;
		}
	} 
	disp(image5,dir,1);
	}
	scanf("%c",name);	

	close_x();

	exit(0);
}	

void disp (unsigned char *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*X_SIZE*3;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/pk%05d.jpg",fram);
	if (ab){jayit(image2,X_SIZE, Y_SIZE, input);}
	free (input);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};


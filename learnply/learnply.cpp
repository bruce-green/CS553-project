/*

Functions for learnply

Eugene Zhang, 2005
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <string.h>
#include <fstream>
#include <vector>
#include "ppm.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "tmatrix.h"
#include "trackball.h"

using namespace std;

#define	NPN 128
#define NMESH  100
#define DM  ((float) (1.0/(NMESH-1.0)))
#define NPIX  512
#define SCALE 4.0
#define M_PI 3.14159265358979323846

int    iframe = 0;
int    Npat = 32;
int    alpha = (0.12 * 255);
float  sa;
float  tmax = NPIX / (SCALE*NPN);
float  dmax = SCALE / NPIX;

int index = 0;

double dsep = 0.1;
double dtest = 0.5 * dsep;

GLdouble wx = 0.5, wy = 0.5, wz;

static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width = 1024;
const int win_height = 1024;

double radius_factor = 0.9;

int display_mode = 0;
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;
int ACSIZE = 1; // for antialiasing
int view_mode = 0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

double zoom = 1.0;
double translation[2] = { 0, 0 };

double xi = 0.4;
double yi = 0.4;

struct Streamline
{
	std::vector<Point>streamline;
};

std::vector<Point>points;
std::vector<Point> samplePoint1;
std::vector<Point> samplePoint2;
std::vector<Point> samplePoint3;
std::vector<Point> samplePoint4;
std::vector<Streamline> samplePoint;
std::vector<Streamline> stlines;
std::vector<Streamline> stlines_b;

Polyhedron *poly;

void init(void);
void initGL(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);
void makePatterns(void);
void display_ibfv(void);
void Advect(double x, double y);
void Advect_b(double x, double y);
void getDP(double x, double y, double *px, double *py);
void seed_point(void);
bool istaken(double x, double y);
bool istaken1(double x, double y);

/******************************************************************************
Main program.
******************************************************************************/
int main(int argc, char** argv)
{
	char *progname;
	int num = 1;
	FILE *this_file;

	points.push_back(Point(xi, yi));

	progname = argv[0];

	this_file = fopen("../new_vector_data/source.ply", "r"); //load ply file
	poly = new Polyhedron(this_file);
	fclose(this_file);
	mat_ident(rotmat);

	poly->initialize(); // initialize everything
	poly->calc_bounding_sphere(); //calculate the coordinate
	poly->calc_face_normals_and_area(); //calculate the normals for the face

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	//glutInitWindowSize(win_width, win_height);
	glutInitWindowSize(NPIX, NPIX);
	//glutCreateWindow(argv[0]);
	glutCreateWindow("Geometric Modeling");
	initGL();

	//makePatterns();

	glutDisplayFunc(display_ibfv);
	glutIdleFunc(display_ibfv);

	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMainLoop();

	poly->finalize();  // finalize everything
	return 0;
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/
Polyhedron::Polyhedron(FILE *file)
{
	int i, j;
	int elem_count;
	char *elem_name;
	/*** Read in the original PLY object ***/
	in_ply = read_ply(file);

	for (i = 0; i < in_ply->num_elem_types; i++) {

		/* prepare to read the i'th list of elements */
		elem_name = setup_element_read_ply(in_ply, i, &elem_count);

		if (equal_strings("vertex", elem_name)) {

			/* create a vertex list to hold all the vertices */
			nverts = max_verts = elem_count;
			vlist = new Vertex *[nverts];

			/* set up for getting vertex elements */

			setup_property_ply(in_ply, &vert_props[0]); //x
			setup_property_ply(in_ply, &vert_props[1]); //y
			setup_property_ply(in_ply, &vert_props[2]); //z
			setup_property_ply(in_ply, &vert_props[3]); //vx
			setup_property_ply(in_ply, &vert_props[4]); //vy
			setup_property_ply(in_ply, &vert_props[5]); //vz
			setup_property_ply(in_ply, &vert_props[6]); //s

			vert_other = get_other_properties_ply(in_ply,
				offsetof(Vertex_io, other_props));

			/* grab all the vertex elements */
			for (j = 0; j < nverts; j++) {
				Vertex_io vert;
				get_element_ply(in_ply, (void *)&vert);

				/* copy info from the "vert" structure */
				vlist[j] = new Vertex(vert.x, vert.y, vert.z);
				vlist[j]->vx = vert.vx;
				vlist[j]->vy = vert.vy;
				vlist[j]->vz = vert.vz;
				vlist[j]->scalar = vert.s;
				vlist[j]->other_props = vert.other_props;
			}
		}
		else if (equal_strings("face", elem_name)) {

			/* create a list to hold all the face elements */
			nquads = max_quads = elem_count;
			qlist = new Quad *[nquads];

			/* set up for getting face elements */
			setup_property_ply(in_ply, &face_props[0]);
			face_other = get_other_properties_ply(in_ply, offsetof(Face_io, other_props));

			/* grab all the face elements */
			for (j = 0; j < elem_count; j++) {
				Face_io face;
				get_element_ply(in_ply, (void *)&face);

				if (face.nverts != 4) {
					fprintf(stderr, "Face has %d vertices (should be four).\n",
						face.nverts);
					exit(-1);
				}

				/* copy info from the "face" structure */
				qlist[j] = new Quad;
				qlist[j]->nverts = 4;
				qlist[j]->verts[0] = (Vertex *)face.verts[0];
				qlist[j]->verts[1] = (Vertex *)face.verts[1];
				qlist[j]->verts[2] = (Vertex *)face.verts[2];
				qlist[j]->verts[3] = (Vertex *)face.verts[3];
				qlist[j]->other_props = face.other_props;
			}
		}
		else
			get_other_element_ply(in_ply);
	}

	/* close the file */
	close_ply(in_ply);

	/* fix up vertex pointers in quads */
	for (i = 0; i < nquads; i++) {
		qlist[i]->verts[0] = vlist[(int)qlist[i]->verts[0]];
		qlist[i]->verts[1] = vlist[(int)qlist[i]->verts[1]];
		qlist[i]->verts[2] = vlist[(int)qlist[i]->verts[2]];
		qlist[i]->verts[3] = vlist[(int)qlist[i]->verts[3]];
	}

	/* get rid of quads that use the same vertex more than once */

	for (i = nquads - 1; i >= 0; i--) {

		Quad *quad = qlist[i];
		Vertex *v0 = quad->verts[0];
		Vertex *v1 = quad->verts[1];
		Vertex *v2 = quad->verts[2];
		Vertex *v3 = quad->verts[3];

		if (v0 == v1 || v1 == v2 || v2 == v3 || v3 == v0) {
			free(qlist[i]);
			nquads--;
			qlist[i] = qlist[nquads];
		}
	}

	//printf("success");
}

/******************************************************************************
Write out a polyhedron to a file.   (no change)
******************************************************************************/
void Polyhedron::write_file(FILE *file)
{
	int i;
	PlyFile *ply;
	char **elist;
	int num_elem_types;

	/*** Write out the transformed PLY object ***/

	elist = get_element_list_ply(in_ply, &num_elem_types);
	ply = write_ply(file, num_elem_types, elist, in_ply->file_type);

	/* describe what properties go into the vertex elements */

	describe_element_ply(ply, "vertex", nverts);
	describe_property_ply(ply, &vert_props[0]);
	describe_property_ply(ply, &vert_props[1]);
	describe_property_ply(ply, &vert_props[2]);
	//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

	describe_element_ply(ply, "face", nquads);
	describe_property_ply(ply, &face_props[0]);

	//  describe_other_properties_ply (ply, face_other,
	//                                offsetof(Face_io,other_props));

	//  describe_other_elements_ply (ply, in_ply->other_elems);

	copy_comments_ply(ply, in_ply);
	char mm[1024];
	printf(mm, "modified by learnply");
	//  append_comment_ply (ply, "modified by simvizply %f");
	append_comment_ply(ply, mm);
	copy_obj_info_ply(ply, in_ply);

	header_complete_ply(ply);

	/* set up and write the vertex elements */
	put_element_setup_ply(ply, "vertex");
	for (i = 0; i < nverts; i++) {
		Vertex_io vert;

		/* copy info to the "vert" structure */
		vert.x = vlist[i]->x;
		vert.y = vlist[i]->y;
		vert.z = vlist[i]->z;
		vert.other_props = vlist[i]->other_props;

		put_element_ply(ply, (void *)&vert);
	}

	/* index all the vertices */
	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	/* set up and write the face elements */
	put_element_setup_ply(ply, "face");

	Face_io face;
	face.verts = new int[4];

	for (i = 0; i < nquads; i++) {

		/* copy info to the "face" structure */
		face.nverts = 4;
		face.verts[0] = qlist[i]->verts[0]->index;
		face.verts[1] = qlist[i]->verts[1]->index;
		face.verts[2] = qlist[i]->verts[2]->index;
		face.verts[3] = qlist[i]->verts[3]->index;
		face.other_props = qlist[i]->other_props;

		put_element_ply(ply, (void *)&face);
	}
	put_other_elements_ply(ply);

	close_ply(ply);
	free_ply(ply);
}

/******************************************************************************
Initialization. (no change)
******************************************************************************/
void Polyhedron::initialize() {
	icVector3 v1, v2; //v1,v2 - two vertices of f1 that define edge

	create_pointers();
	calc_edge_length();
	seed = -1;
}

/******************************************************************************
Finalization. (no change)
******************************************************************************/
void Polyhedron::finalize() {
	int i;

	for (i = 0; i < nquads; i++) {
		free(qlist[i]->other_props);
		free(qlist[i]);
	}
	for (i = 0; i < nedges; i++) {
		free(elist[i]->quads);
		free(elist[i]);
	}
	for (i = 0; i < nverts; i++) {
		free(vlist[i]->quads);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	free(qlist); //clear up the list
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face      (no change)
******************************************************************************/
Quad *Polyhedron::find_common_edge(Quad *f1, Vertex *v1, Vertex *v2)
{
	int i, j;
	Quad *f2;
	Quad *adjacent = NULL;

	/* look through all faces of the first vertex */

	for (i = 0; i < v1->nquads; i++) {
		f2 = v1->quads[i];
		if (f2 == f1)
			continue;
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < f2->nverts; j++) {

			/* look for a match */
			if (f2->verts[j] == v2) {

#if 0
				/* watch out for triple edges */

				if (adjacent != NULL) {

					fprintf(stderr, "model has triple edges\n");

					fprintf(stderr, "face 1: ");
					for (k = 0; k < f1->nverts; k++)
						fprintf(stderr, "%d ", f1->iverts[k]);
					fprintf(stderr, "\nface 2: ");
					for (k = 0; k < f2->nverts; k++)
						fprintf(stderr, "%d ", f2->iverts[k]);
					fprintf(stderr, "\nface 3: ");
					for (k = 0; k < adjacent->nverts; k++)
						fprintf(stderr, "%d ", adjacent->iverts[k]);
					fprintf(stderr, "\n");

				}

				/* if we've got a match, remember this face */
				adjacent = f2;
#endif

#if 1
				/* if we've got a match, return this face */
				return (f2);
#endif

			}
		}
	}

	return (adjacent);
}

/******************************************************************************
Create an edge. (no change)

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/
void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
	int i, j;
	Quad *f;

	/* make sure there is enough room for a new edge */

	if (nedges >= max_edges) {

		max_edges += 100;
		Edge **list = new Edge *[max_edges];

		/* copy the old list to the new one */
		for (i = 0; i < nedges; i++)
			list[i] = elist[i];

		/* replace list */
		free(elist);
		elist = list;
	}

	/* create the edge */

	elist[nedges] = new Edge;
	Edge *e = elist[nedges];
	e->index = nedges;
	e->verts[0] = v1;
	e->verts[1] = v2;
	nedges++;

	/* count all quads that will share the edge, and do this */
	/* by looking through all faces of the first vertex */

	e->nquads = 0;

	for (i = 0; i < v1->nquads; i++) {
		f = v1->quads[i];
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 4; j++) {
			/* look for a match */
			if (f->verts[j] == v2) {
				e->nquads++;
				break;
			}
		}
	}

	/* make room for the face pointers (at least two) */
	if (e->nquads < 2)
		e->quads = new Quad *[2];
	else
		e->quads = new Quad *[e->nquads];

	/* create pointers from edges to faces and vice-versa */

	e->nquads = 0; /* start this out at zero again for creating ptrs to quads */

	for (i = 0; i < v1->nquads; i++) {

		f = v1->quads[i];

		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 4; j++)
			if (f->verts[j] == v2) {

				e->quads[e->nquads] = f;
				e->nquads++;

				if (f->verts[(j + 1) % 4] == v1)
					f->edges[j] = e;
				else if (f->verts[(j + 2) % 4] == v1)
					f->edges[(j + 2) % 4] = e;
				else if (f->verts[(j + 3) % 4] == v1)
					f->edges[(j + 3) % 4] = e;
				else {
					fprintf(stderr, "Non-recoverable inconsistancy in create_edge()\n");
					exit(-1);
				}

				break;  /* we'll only find one instance of v2 */
			}

	}
}

/******************************************************************************
Create edges (no change)
******************************************************************************/
void Polyhedron::create_edges()
{
	int i, j;
	Quad *f;
	Vertex *v1, *v2;
	double count = 0;

	/* count up how many edges we may require */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < f->nverts; j++) {
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			Quad *result = find_common_edge(f, v1, v2);
			if (result)
				count += 0.5;
			else
				count += 1;
		}
	}

	/*
	printf ("counted %f edges\n", count);
	*/

	/* create space for edge list */

	max_edges = (int)(count + 10);  /* leave some room for expansion */
	elist = new Edge *[max_edges];
	nedges = 0;

	/* zero out all the pointers from faces to edges */

	for (i = 0; i < nquads; i++)
		for (j = 0; j < 4; j++)
			qlist[i]->edges[j] = NULL;

	/* create all the edges by examining all the quads */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < 4; j++) {
			/* skip over edges that we've already created */
			if (f->edges[j])
				continue;
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			create_edge(v1, v2);
		}
	}
}

/******************************************************************************
Create pointers from vertices to faces. (no change)
******************************************************************************/
void Polyhedron::vertex_to_quad_ptrs()
{
	int i, j;
	Quad *f;
	Vertex *v;

	/* zero the count of number of pointers to faces */

	for (i = 0; i < nverts; i++)
		vlist[i]->max_quads = 0;

	/* first just count all the face pointers needed for each vertex */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < f->nverts; j++)
			f->verts[j]->max_quads++;
	}

	/* allocate memory for face pointers of vertices */

	for (i = 0; i < nverts; i++) {
		vlist[i]->quads = (Quad **)
			malloc(sizeof(Quad *) * vlist[i]->max_quads);
		vlist[i]->nquads = 0;
	}

	/* now actually create the face pointers */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < f->nverts; j++) {
			v = f->verts[j];
			v->quads[v->nquads] = f;
			v->nquads++;
		}
	}
}

/******************************************************************************
Find the other quad that is incident on an edge, or NULL if there is
no other. (no change)
******************************************************************************/
Quad *Polyhedron::other_quad(Edge *edge, Quad *quad)
{
	/* search for any other quad */
	for (int i = 0; i < edge->nquads; i++)
		if (edge->quads[i] != quad)
			return (edge->quads[i]);

	/* there is no such other quad if we get here */
	return (NULL);
}

/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered  (no change)
******************************************************************************/
void Polyhedron::order_vertex_to_quad_ptrs(Vertex *v)
{
	int i, j;
	Quad *f;
	Quad *fnext;
	int nf;
	int vindex;
	int boundary;
	int count;

	nf = v->nquads;
	f = v->quads[0];

	/* go backwards (clockwise) around faces that surround a vertex */
	/* to find out if we reach a boundary */

	boundary = 0;

	for (i = 1; i <= nf; i++) {

		/* find reference to v in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[j] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #1\n");
			exit(-1);
		}

		/* corresponding face is the previous one around v */
		fnext = other_quad(f->edges[vindex], f);

		/* see if we've reached a boundary, and if so then place the */
		/* current face in the first position of the vertice's face list */

		if (fnext == NULL) {
			/* find reference to f in v */
			for (j = 0; j < v->nquads; j++)
				if (v->quads[j] == f) {
					v->quads[j] = v->quads[0];
					v->quads[0] = f;
					break;
				}
			boundary = 1;
			break;
		}

		f = fnext;
	}

	/* now walk around the faces in the forward direction and place */
	/* them in order */

	f = v->quads[0];
	count = 0;

	for (i = 1; i < nf; i++) {

		/* find reference to vertex in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[(j + 1) % f->nverts] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #2\n");
			exit(-1);
		}

		/* corresponding face is next one around v */
		fnext = other_quad(f->edges[vindex], f);

		/* break out of loop if we've reached a boundary */
		count = i;
		if (fnext == NULL) {
			break;
		}

		/* swap the next face into its proper place in the face list */
		for (j = 0; j < v->nquads; j++)
			if (v->quads[j] == fnext) {
				v->quads[j] = v->quads[i];
				v->quads[i] = fnext;
				break;
			}

		f = fnext;
	}
}

/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found  (no change)
******************************************************************************/
int Polyhedron::face_to_vertex_ref(Quad *f, Vertex *v)
{
	int j;
	int vindex = -1;

	for (j = 0; j < f->nverts; j++)
		if (f->verts[j] == v) {
			vindex = j;
			break;
		}

	return (vindex);
}

/******************************************************************************
Create various face and vertex pointers. (no change)
******************************************************************************/
void Polyhedron::create_pointers()
{
	int i;

	/* index the vertices and quads */

	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	for (i = 0; i < nquads; i++)
		qlist[i]->index = i;

	/* create pointers from vertices to quads */
	vertex_to_quad_ptrs();

	/* make edges */
	create_edges();


	/* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
		order_vertex_to_quad_ptrs(vlist[i]);

	}
	/* index the edges */

	for (i = 0; i < nedges; i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
		elist[i]->index = i;
	}

}

/******************************************************************************
For this cs453/553, this function is used to calculate the coordinates,since the
original cooridnates in ply file does not equal to the coordinate in the world
coordinate in this program (no change)
******************************************************************************/
void Polyhedron::calc_bounding_sphere()
{
	unsigned int i;
	icVector3 min, max;

	for (i = 0; i < nverts; i++) {
		if (i == 0) {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		}
		else {
			if (vlist[i]->x < min.entry[0])
				min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
				max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
				min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
				max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
				min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
				max.entry[2] = vlist[i]->z;
		}
	}
	center = (min + max) * 0.5;
	radius = length(center - min);
}

/******************************************************************************
Calculate the length of the edge (no change)
******************************************************************************/
void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i = 0; i < nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1 - v2);
	}
}

/******************************************************************************
Calculate the normals and area (no change)
******************************************************************************/
void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2, v3;
	Quad *temp_q;
	double edge_length[4];

	area = 0.0;
	for (i = 0; i < nquads; i++) {
		for (j = 0; j < 4; j++)
			edge_length[j] = qlist[i]->edges[j]->length;

		icVector3 d1, d2;
		d1.set(qlist[i]->verts[0]->x, qlist[i]->verts[0]->y, qlist[i]->verts[0]->z);
		d2.set(qlist[i]->verts[2]->x, qlist[i]->verts[2]->y, qlist[i]->verts[2]->z);
		double dia_length = length(d1 - d2);

		double temp_s1 = (edge_length[0] + edge_length[1] + dia_length) / 2.0;
		double temp_s2 = (edge_length[2] + edge_length[3] + dia_length) / 2.0;
		qlist[i]->area = sqrt(temp_s1*(temp_s1 - edge_length[0])*(temp_s1 - edge_length[1])*(temp_s1 - dia_length)) +
			sqrt(temp_s2*(temp_s2 - edge_length[2])*(temp_s2 - edge_length[3])*(temp_s2 - dia_length));

		area += qlist[i]->area;
		temp_q = qlist[i];
		v1.set(vlist[qlist[i]->verts[0]->index]->x, vlist[qlist[i]->verts[0]->index]->y, vlist[qlist[i]->verts[0]->index]->z);
		v2.set(vlist[qlist[i]->verts[1]->index]->x, vlist[qlist[i]->verts[1]->index]->y, vlist[qlist[i]->verts[1]->index]->z);
		v0.set(vlist[qlist[i]->verts[2]->index]->x, vlist[qlist[i]->verts[2]->index]->y, vlist[qlist[i]->verts[2]->index]->z);
		qlist[i]->normal = cross(v0 - v1, v2 - v1);
		normalize(qlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i = 0; i < nquads; i++) {
		icVector3 cent(vlist[qlist[i]->verts[0]->index]->x, vlist[qlist[i]->verts[0]->index]->y, vlist[qlist[i]->verts[0]->index]->z);
		signedvolume += dot(test - cent, qlist[i]->normal)*qlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume < 0)
		orientation = 0;
	else {
		orientation = 1;
		for (i = 0; i < nquads; i++)
			qlist[i]->normal *= -1.0;
	}
}

/******************************************************************************
OpenGL set up and initialization (no change)
******************************************************************************/
void init(void) {
	/* select clearing color */

	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

void initGL(void) {
	// select clearing color 

	glClearColor(0.0, 0.0, 0.0, 0.0);  // background 
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);

	glViewport(0, 0, (GLsizei)NPIX, (GLsizei)NPIX);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '0':
		display_mode = 0;
		display();
		break;

	case '1':
		display_mode = 1;
		display();
		break;

	case '2':
		display_mode = 2;
		display();
		break;

	case '3':
		display_mode = 3;
		display();
		break;

	case '4':
		display_mode = 4;
		display();
		break;

	case '5':
		display_mode = 5;
		display();
		break;

	case '6':
		display_mode = 6;
		display();
		break;

	case '7':
		display_mode = 7;
		display();
		break;

	case '8':
		display_mode = 8;
		display();
		break;

	case '9':
		display_mode = 9;
		display();
		break;

	case '|':
		this_file = fopen("rotmat.txt", "w");
		for (i = 0; i < 4; i++)
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
		break;

	case '^':
		this_file = fopen("rotmat.txt", "r");
		for (i = 0; i < 4; i++)
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
		break;

	case 'r':
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		display();
		break;

	}
}

/******************************************************************************
define the number of vertices (no change)
******************************************************************************/
Polyhedron::Polyhedron()
{
	nverts = nedges = nquads = 0;
	max_verts = max_quads = 50;

	vlist = new Vertex *[max_verts];
	qlist = new Quad *[max_quads];
}

/******************************************************************************
openGL setup (no change)
******************************************************************************/
void multmatrix(const Matrix m)
{
	int i, j, index = 0;

	GLfloat mat[16];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf(mat);
}

/******************************************************************************
openGL view setup (light) (no change)
******************************************************************************/
void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor * zoom, radius_factor*zoom, -radius_factor * zoom, radius_factor*zoom, 0.0, 40.0);
	else
		glFrustum(-radius_factor * zoom, radius_factor*zoom, -radius_factor, radius_factor, -1000, 1000);
	//gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
openGL view setup (position) (no change)
******************************************************************************/
void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(translation[0], translation[1], -3.0);
	multmatrix(rotmat);

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Assign a name for each element (no change)
******************************************************************************/
int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	//printf("hits = %d\n", hits);
	ptr = (GLuint *)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	//printf("Quad id = %d\n", seed_id);
	return seed_id;
}

/******************************************************************************
Mouse function (no change)
******************************************************************************/
void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;


	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case -1:

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 3:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;

	case 2:

		zoom *= exp(2.0*(t - t_old));

		s_old = s;
		t_old = t;

		display();
		break;

	}
}

/******************************************************************************
Mouse function (no change)
******************************************************************************/
void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON) {
		switch (mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float)win_width;
				float ysize = (float)win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_button = button;
				last_x = x;
				last_y = y;

				/*rotate*/
				if (button == GLUT_RIGHT_BUTTON)
				{
					mouse_mode = -1;
				}

				/*translate*/
				if (button == GLUT_LEFT_BUTTON)
				{
					mouse_mode = 3;
				}


				if (key == GLUT_ACTIVE_SHIFT) {

					GLint viewport[4];
					GLdouble mvmatrix[16], projmatrix[16];
					GLint realy;

					Points temp_p;

					temp_p.index = 1;

					glGetIntegerv(GL_VIEWPORT, viewport);
					glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
					realy = viewport[3] - y - 1;
					gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x,
						&temp_p.y, &temp_p.z);
					//points.push_back(temp_p);

				}
				if (key == GLUT_ACTIVE_ALT) {

					GLint viewport[4];
					GLdouble mvmatrix[16], projmatrix[16];
					GLint realy;

					Points temp_p;

					temp_p.index = 2;

					glGetIntegerv(GL_VIEWPORT, viewport);
					glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
					realy = viewport[3] - y - 1;
					gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x,
						&temp_p.y, &temp_p.z);

					//points.push_back(temp_p);

				}
				if (key == GLUT_ACTIVE_CTRL) {

					GLint viewport[4];
					GLdouble mvmatrix[16], projmatrix[16];
					GLint realy;

					Points temp_p;

					temp_p.index = 3;

					glGetIntegerv(GL_VIEWPORT, viewport);
					glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
					realy = viewport[3] - y - 1;
					gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x,
						&temp_p.y, &temp_p.z);

					//points.push_back(temp_p);

				}
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	}
	else if (button == GLUT_RIGHT_BUTTON) {
		switch (mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float)win_width;
				float ysize = (float)win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_button = button;
				last_x = x;
				last_y = y;

				/*rotate*/
				if (button == GLUT_RIGHT_BUTTON)
				{
					mouse_mode = -1;
				}

				/*translate*/
				if (button == GLUT_LEFT_BUTTON)
				{
					mouse_mode = 3;
				}


				if (key == GLUT_ACTIVE_SHIFT) {

					GLint viewport[4];
					GLdouble mvmatrix[16], projmatrix[16];
					GLint realy;

					Points temp_p;

					temp_p.index = 4;

					glGetIntegerv(GL_VIEWPORT, viewport);
					glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
					realy = viewport[3] - y - 1;
					gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x,
						&temp_p.y, &temp_p.z);

					//points.push_back(temp_p);

				}
				if (key == GLUT_ACTIVE_ALT) {

					GLint viewport[4];
					GLdouble mvmatrix[16], projmatrix[16];
					GLint realy;

					Points temp_p;

					temp_p.index = 5;

					glGetIntegerv(GL_VIEWPORT, viewport);
					glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
					realy = viewport[3] - y - 1;
					//printf("Coordinate at curosr are (%4d, %4d, %4d)\n", x, y, viewport[3]);
					gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x,
						&temp_p.y, &temp_p.z);
					//printf("World coords at z=0 are (%f, %f, %f)\n", temp_p.x, temp_p.y, temp_p.z);
					//points.push_back(temp_p);
					//gluUnProject(x, realy, 1, mvmatrix, projmatrix, viewport, &wx, &wy, &wz);
					//printf("World coords at z=1 are (%f, %f, %f)\n", wx, wy, wz);
				}
				if (key == GLUT_ACTIVE_CTRL) {

					GLint viewport[4];
					GLdouble mvmatrix[16], projmatrix[16];
					GLint realy;

					Points temp_p;

					temp_p.index = 6;

					glGetIntegerv(GL_VIEWPORT, viewport);
					glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
					glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
					realy = viewport[3] - y - 1;
					//printf("Coordinate at curosr are (%4d, %4d, %4d)\n", x, y, viewport[3]);
					gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x,
						&temp_p.y, &temp_p.z);
					//printf("World coords at z=0 are (%f, %f, %f)\n", temp_p.x, temp_p.y, temp_p.z);
					//points.push_back(temp_p);
					//gluUnProject(x, realy, 1, mvmatrix, projmatrix, viewport, &wx, &wy, &wz);
					//printf("World coords at z=1 are (%f, %f, %f)\n", wx, wy, wz);
				}
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	}
	else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLint viewport[4];
			GLdouble mvmatrix[16], projmatrix[16];
			GLint realy;

			Points temp_p;

			glGetIntegerv(GL_VIEWPORT, viewport);
			glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
			glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
			realy = viewport[3] - y - 1;
			gluUnProject(x, realy, 0, mvmatrix, projmatrix, viewport, &temp_p.x, &temp_p.y, &temp_p.z);
			//points.push_back(temp_p);
		}
	}
}

/******************************************************************************
Different visualizations of case
******************************************************************************/
void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	float colorX, colorY;
	double L = (poly->radius * 2) / 30;

	//openGL setting including the materials, lighting
	for (i = 0; i < this_poly->nquads; i++) { //go through all the quads
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad *temp_q = this_poly->qlist[i];

		switch (display_mode) {
		case 0:
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glEnable(GL_LIGHT1);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			}
			else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex *temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i == this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 1:    //draw mesh 
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //"GL_LINE" is for mesh drawing
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex *temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 2:   //checker board 
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex *temp_v = temp_q->verts[j];

				colorX = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				colorY = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				glColor3f(colorX, colorY, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 3:  //grey scale (fx)
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex *temp_v = temp_q->verts[j];
				double new_color = temp_v->f_x;

				glColor3f(new_color, new_color, new_color);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 4:  //bi_color scale (fx)
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);

			for (j = 0; j < 4; j++) {
				Vertex *temp_v = temp_q->verts[j];
				icVector3 c1(1.0f, 0.0f, 0.0f); //red
				icVector3 c2(0.0f, 0.0f, 1.0f); //blue
				icVector3 c = c1 * temp_v->f_x + c2 * temp_v->f_x1;
				glColor3f(c.x, c.y, c.z);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 5:  //height with grey scale (fx)
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);

			for (j = 0; j < 4; j++) {
				Vertex *temp_v = temp_q->verts[j];
				double new_color = temp_v->f_x;
				glColor3f(new_color, new_color, new_color);
				glVertex3d(temp_v->x, temp_v->y, 50 * temp_v->f_x);
			}
			glEnd();
			break;


		case 6:  //grey scale (gx), if you need g(x) for other question, just replace the f_x with g_x
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex *temp_v = temp_q->verts[j];
				double new_color = temp_v->g_x;

				glColor3f(new_color, new_color, new_color);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}

	}
}

/******************************************************************************
Display function
******************************************************************************/
void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv(GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix();
		set_scene(GL_RENDER, poly);
		display_shape(GL_RENDER, poly);
		glPopMatrix();
		glAccum(GL_ACCUM, 1.0 / ACSIZE);
	}
	glAccum(GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}


unsigned char *pixels;
void display_ibfv(void)
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_DEPTH_TEST);

	glTexParameteri(GL_TEXTURE_2D,
		GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,
		GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,
		GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,
		GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV,
		GL_TEXTURE_ENV_MODE, GL_REPLACE);


	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	int   i, j, k;
	double x1, x2, y, px, py;

	for (int k = 0; k < points.size(); k++)
	{
		Point temp_ps = points[k];
		Advect(temp_ps.x, temp_ps.y);

		glLineWidth(2.0);

		for (int i = 0; i < stlines.size(); i++)
		{
			Streamline temp_stline = stlines[i];
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j < stlines[i].streamline.size() - 1; j++) {
				Point temp_points1 = stlines[i].streamline[j];
				Point temp_points2 = stlines[i].streamline[j + 1];
				if (temp_points1.x >= 1 || temp_points1.x <= 0 || temp_points2.x >= 1 || temp_points2.x <= 0) {
					break;
				}
				if (temp_points1.y >= 1 || temp_points1.y <= 0 || temp_points2.y >= 1 || temp_points2.y <= 0) {
					break;
				}
				glColor3f(1.0, 0.0, 0.0);
				glVertex2d(temp_points2.x, temp_points2.y);
			}
			glEnd();
		}
		stlines.clear();

	}

	for (int k = 0; k < points.size(); k++)
	{
		Point temp_ps = points[k];
		Advect_b(temp_ps.x, temp_ps.y);

		glLineWidth(2.0);

		for (int i = 0; i < stlines_b.size(); i++)
		{
			Streamline temp_stline = stlines_b[i];
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j < stlines_b[i].streamline.size() - 1; j++) {
				Point temp_points1 = stlines_b[i].streamline[j];
				Point temp_points2 = stlines_b[i].streamline[j + 1];
				//printf("x:%f, y:%f\n", temp_points.x, temp_points.y);
				if (temp_points1.x >= 1 || temp_points1.x <= 0 || temp_points2.x >= 1 || temp_points2.x <= 0) {
					break;
				}
				if (temp_points1.y >= 1 || temp_points1.y <= 0 || temp_points2.y >= 1 || temp_points2.y <= 0) {
					break;
				}
				glColor3f(1.0, 0.0, 0.0);
				glVertex2d(temp_points2.x, temp_points2.y);
			}
			glEnd();
		}
		stlines_b.clear();
	}

	seed_point();

	if (index < samplePoint.size()) {
		index += 2;
	}

	glutSwapBuffers();
}

void makePatterns(void)
{
	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k, t;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++) {
		for (j = 0; j < NPN; j++) {
			phase[i][j] = rand() % 256;
		}
	}

	for (k = 0; k < Npat; k++) {
		t = k * 256 / Npat;
		for (i = 0; i < NPN; i++) {
			for (j = 0; j < NPN; j++) {
				/*pat[i][j][0] = ppms.r[(NPN - i - 1) * NPN + j];
				pat[i][j][1] = ppms.g[(NPN - i - 1) * NPN + j];
				pat[i][j][2] = ppms.b[(NPN - i - 1) * NPN + j];
				pat[i][j][3] = alpha;*/
				pat[i][j][0] =
					pat[i][j][1] =
					//pat[i][j][2] = lut[(t + phase[i][j]) % 255];
					pat[i][j][2] = 0;
				pat[i][j][3] = alpha;
			}
			glNewList(k + 1, GL_COMPILE);
			glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
			glEndList();
		}
	}
}


void getDP(double x, double y, double *px, double *py) {
	double dx, dy, vx, vy, rx, ry, r;
	dx = -(x - 0.5);
	dy = (y - 0.5);

	r = dx * dx + dy * dy;

	if (r > dmax*dmax) {
		r = sqrt(r);
		dx *= dmax / r;
		dy *= dmax / r;
	}
	*px = dx;
	*py = dy;
}

void getDP_b(double x, double y, double *px, double *py) {
	double dx, dy, vx, vy, rx, ry, r;
	dx = (x - 0.5);
	dy = -(y - 0.5);

	r = dx * dx + dy * dy;

	if (r > dmax*dmax) {
		r = sqrt(r);
		dx *= dmax / r;
		dy *= dmax / r;
	}
	*px = dx;
	*py = dy;
}


void Advect(double x, double y)
{
	double step = 0.6;
	double x0 = x;
	double y0 = y;

	Streamline streamlines1;
	Streamline samplePoints;
	double d0,d1;
	
	streamlines1.streamline.push_back(Point(x0, y0));
	samplePoints.streamline.push_back(Point(x0, y0));
	if (x0 > 0.5 && y0 > 0.5) {
		samplePoint1.push_back(Point(x0, y0));
	}
	else if (x0 < 0.5 && y0 > 0.5)
	{
		samplePoint2.push_back(Point(x0, y0));
	}
	else if (x0 < 0.5 && y0 < 0.5)
	{
		samplePoint3.push_back(Point(x0, y0));
	}
	else if (x0 > 0.5 && y0 < 0.5)
	{
		samplePoint4.push_back(Point(x0, y0));
	}
	
	for (int i = 0; i < 100; i++) {
		double dx0, dy0;
		getDP(x0, y0, &dx0, &dy0);
		double x1 = x0 + step * dx0;
		double y1 = y0 + step * dy0;

		double dx1, dy1;
		getDP(x1, y1, &dx1, &dy1);
		double vx = (dx0 + dx1) / 2;
		double vy = (dy0 + dy1) / 2;

		double x2 = x0 + step * vx;
		double y2 = y0 + step * vy;

		streamlines1.streamline.push_back(Point(x2, y2));
		if (i % 5 == 0)
		{
			samplePoints.streamline.push_back(Point(x2, y2));
			if (x0 > 0.5 && y0 > 0.5) {
				samplePoint1.push_back(Point(x2, y2));
			}
			else if (x0 < 0.5 && y0 > 0.5)
			{
				samplePoint2.push_back(Point(x2, y2));
			}
			else if (x0 < 0.5 && y0 < 0.5)
			{
				samplePoint3.push_back(Point(x2, y2));
			}
			else if (x0 > 0.5 && y0 < 0.5)
			{
				samplePoint4.push_back(Point(x2, y2));
			}
		}
	
		x0 = x2;
		y0 = y2;
		
	}
	stlines.push_back(streamlines1);
	samplePoint.push_back(samplePoints);
}

void Advect_b(double x, double y)
{
	double step = 0.6;
	double x0 = x;
	double y0 = y;

	Streamline streamlines2;
	Streamline samplePoints;
		
	streamlines2.streamline.push_back(Point(x0, y0));
	samplePoints.streamline.push_back(Point(x0, y0));
	if (x0 > 0.5 && y0 > 0.5) {
		samplePoint1.push_back(Point(x0, y0));
	}
	else if (x0 < 0.5 && y0 > 0.5)
	{
		samplePoint2.push_back(Point(x0, y0));
	}
	else if (x0 < 0.5 && y0 < 0.5)
	{
		samplePoint3.push_back(Point(x0, y0));
	}
	else if (x0 > 0.5 && y0 < 0.5)
	{
		samplePoint4.push_back(Point(x0, y0));
	}

	for (int i = 0; i < 100; i++) {
		double dx0, dy0;
		getDP_b(x0, y0, &dx0, &dy0);
		double x1 = x0 + step * dx0;
		double y1 = y0 + step * dy0;

		double dx1, dy1;
		getDP_b(x1, y1, &dx1, &dy1);
		double vx = (dx0 + dx1) / 2;
		double vy = (dy0 + dy1) / 2;

		double x2 = x0 + step * vx;
		double y2 = y0 + step * vy;

		streamlines2.streamline.push_back(Point(x2, y2));
		if (i % 5 == 0)
		{
			samplePoints.streamline.push_back(Point(x2, y2));
			if (x0 > 0.5 && y0 > 0.5) {
				samplePoint1.push_back(Point(x2, y2));
			}
			else if (x0 < 0.5 && y0 > 0.5)
			{
				samplePoint2.push_back(Point(x2, y2));
			}
			else if (x0 < 0.5 && y0 < 0.5)
			{
				samplePoint3.push_back(Point(x2, y2));
			}
			else if (x0 > 0.5 && y0 < 0.5)
			{
				samplePoint4.push_back(Point(x2, y2));
			}
		}

		x0 = x2;
		y0 = y2;
	}
	stlines_b.push_back(streamlines2);
	samplePoint.push_back(samplePoints);
}

void seed_point(void) 
{
	double x, y;
	double d;
	int count = 0;

	Streamline temp_samplePoint;

	if (samplePoint.size() > 0)
	{
		temp_samplePoint = samplePoint[index];
		for (int i = 0; i < temp_samplePoint.streamline.size(); i++)
		{
			Point sPoint = temp_samplePoint.streamline[i];

			for (int j = ((sPoint.x - (2 * dsep)) * NMESH); j < ((sPoint.x + (2 * dsep)) * NMESH); j++)
			{
				x = DM * j;

				for (int k = ((sPoint.y - (2 * dsep)) * NMESH); k < ((sPoint.y + (2 * dsep)) * NMESH); k++)
				{
					y = DM * k;
					d = sqrt((x - sPoint.x) * (x - sPoint.x) + (y - sPoint.y) * (y - sPoint.y));
					
					if (count < 2)
					{
						if (d < (dsep + 0.0001) && d > (dsep - 0.0001) && istaken(x, y) == 0)
						{
							points.push_back(Point(x, y));
							count++;
						}
					}
				}
			}
		}
		
	}
}


bool istaken1(double x, double y) {

	double d1;
	for (int i = 0; i < samplePoint.size(); i++)
	{
		Streamline temp_sp = samplePoint[i];
		for (int j = 0; j < temp_sp.streamline.size(); j++)
		{
			Point sample_p = temp_sp.streamline[j];
			d1 = sqrt((x - sample_p.x) * (x - sample_p.x) + (y - sample_p.y) * (y - sample_p.y));
			if (d1 < dtest)
			{
				return true;
			}
		}
	}
	return false;
}

bool istaken(double x , double y) {
	double d1;
	for (int i = 0; i < points.size(); i++)
	{
		d1 = sqrt((x - points[i].x) * (x - points[i].x) + (y - points[i].y) * (y - points[i].y));
		if (d1 < dtest)
		{
			return true;
		}
	}
	if (x > 0.5 && y > 0.5) {
		for (int i = 0; i < samplePoint1.size(); i++)
		{
			d1 = sqrt((x - samplePoint1[i].x) * (x - samplePoint1[i].x) + (y - samplePoint1[i].y) * (y - samplePoint1[i].y));
			if (d1 < dtest)
			{
				return true;
			}
		}
		return false;
	}
	else if (x < 0.5 && y > 0.5)
	{
		for (int i = 0; i < samplePoint2.size(); i++)
		{
			d1 = sqrt((x - samplePoint2[i].x) * (x - samplePoint2[i].x) + (y - samplePoint2[i].y) * (y - samplePoint2[i].y));
			if (d1 < dtest)
			{
				return true;
			}
		}
		return false;
	}
	else if (x < 0.5 && y < 0.5)
	{
		for (int i = 0; i < samplePoint3.size(); i++)
		{
			d1 = sqrt((x - samplePoint3[i].x) * (x - samplePoint3[i].x) + (y - samplePoint3[i].y) * (y - samplePoint3[i].y));
			if (d1 < dtest)
			{
				return true;
			}
		}
		return false;
	}
	else if (x > 0.5 && y < 0.5)
	{
		for (int i = 0; i < samplePoint4.size(); i++)
		{
			d1 = sqrt((x - samplePoint4[i].x) * (x - samplePoint4[i].x) + (y - samplePoint4[i].y) * (y - samplePoint4[i].y));
			if (d1 < dtest)
			{
				return true;
			}
		}
		return false;
	}
	
}

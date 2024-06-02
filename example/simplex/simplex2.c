// #include <p4est_to_p8est.h>  // This is only for my vim-language-server

#include <sc_options.h>

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <p4est_simplex_mesh.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <p8est_simplex_mesh.h>
#endif

static char *p4est_connectivity_names[] = {
#ifndef P4_TO_P8
  "brick23","corner","cubed","disk","icosahedron","moebius",
  "periodic","pillow","rotwrap","star","shell2d","disk2d",
  "bowtie","unit","sphere2d",
#else
  "brick235","periodic","rotcubes","rotwrap","shell","sphere",
  "twocubes","twowrap","unit",
#endif
  NULL
};

typedef struct
{
  int minlevel, maxlevel;
  const char *conn;
  const char *mshpath;

  int hangine_node_level;
}
options_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  options_t *opts;

  p4est_t *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost_layer;
  p4est_geometry_t *geom;

  char *errmsg;
}
context_t;


#ifdef P4_TO_P8
#define QUAD_COORDS(quad) \
            (quad).x / (double) P4EST_ROOT_LEN, \
            (quad).y / (double) P4EST_ROOT_LEN, \
            (quad).z / (double) P4EST_ROOT_LEN

#else
#define QUAD_COORDS(quad) \
            (quad).x / (double) P4EST_ROOT_LEN, \
            (quad).y / (double) P4EST_ROOT_LEN
#endif

/** Refines the p4est, producing hanging-nodes
 */
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  context_t *g = (context_t *) p4est->user_pointer;
  options_t *opts = g->opts;

  int refine_level = g->opts->maxlevel;
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 1;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    // return 0;
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
#ifndef P4_TO_P8
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }
#else
  if (quadrant->z >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }
#endif

  return 1;
}

int
create_context_from_opts(context_t *g, options_t *opts)
{
  p4est_t *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost_layer;
  p4est_geometry_t *geom = NULL;

  /* Create p4est */
#ifndef P4_TO_P8
  if (strcmp(opts->conn, "sphere2d") == 0) {
    conn = p4est_connectivity_new_cubed();
  }
  else
#endif
  {
    conn = p4est_connectivity_new_byname (opts->conn);
  }

  if (!conn) {
      printf("ERROR: HEY\n");
  }

  g->conn = conn;
  p4est = p4est_new_ext (g->mpicomm, conn, 0, opts->minlevel, 1, 0, NULL, NULL);
  g->p4est = p4est;
  p4est->user_pointer = (void *) g;
  p4est_refine_ext (p4est, 1, opts->maxlevel, refine_fn, NULL, NULL);
  p4est_partition (p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  /* Create ghost layer */
  ghost_layer = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
  g->ghost_layer = ghost_layer;

  /* Create geometry */
#ifndef P4_TO_P8
  if (opts->conn && !strcmp(opts->conn, "disk2d")) {
    geom = p4est_geometry_new_disk2d (conn, 1.0, 2.0);
  } else if (opts->conn && !strcmp(opts->conn, "icosahedron")) {
    geom = p4est_geometry_new_icosahedron (conn, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "shell2d")) {
    geom = p4est_geometry_new_shell2d (conn, 2.0, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "sphere2d")) {
    geom = p4est_geometry_new_sphere2d(conn, 1.0);
  }
#else
  if (opts->conn && !strcmp(opts->conn, "shell")) {
    geom = p8est_geometry_new_shell (conn, 2.0, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "sphere")) {
    geom = p8est_geometry_new_sphere (conn, 3.0, 2.0, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "torus")) {
    geom = p8est_geometry_new_torus (conn, 1.0, 2.0, 4.0);
  }
#endif
  else {
    geom = p4est_geometry_new_connectivity(conn);
  }
  g->geom = geom;

  return 0;
}

void
destroy_context(context_t *g)
{
  if (g->p4est)
    p4est_destroy(g->p4est);
  if (g->conn)
    p4est_connectivity_destroy(g->conn);
  if (g->ghost_layer)
    p4est_ghost_destroy(g->ghost_layer);
  if (g->geom)
    p4est_geometry_destroy(g->geom);
}


int
main(int argc, char **argv) {
  int mpiret;
  int first_argc;
  int exitcode;
  sc_options_t *sc_opts;
  options_t options, *opts = &options;
  context_t context = {0}, *g = &context;
  p4est_simplex_mesh_t *smesh;
  char filepath[1024];


  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  g->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (g->mpicomm, &g->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->mpicomm, &g->mpirank);
  SC_CHECK_MPI (mpiret);

  /* Does this really need to be before the parsing? */
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  // User info flags
  int help = 0;
  int list_conn = 0;

  smesh = NULL;

  /* Initialize default option values */
  opts->minlevel = 2;
  opts->maxlevel = 5;
  opts->conn = "unit";
  opts->mshpath = NULL;  // Defaults to "{opts->conn}.msh"

  opts->hangine_node_level = 2;

  sc_opts = sc_options_new (argv[0]);
  sc_options_add_int (sc_opts, 'l', "minlevel", &opts->minlevel, opts->minlevel, "Lowest level");
  sc_options_add_int (sc_opts, 'L', "maxlevel", &opts->maxlevel, opts->maxlevel, "Highest level");

  sc_options_add_int (sc_opts, 'h', "hanging-node-level", &opts->hangine_node_level, opts->hangine_node_level, 
      "This determines what hanging nodes to include. h=0,1,2\n"
      // " h=-1: Do not add element centre nodes either,\n"
      " h=0: Include no hanging nodes,\n"
      " h=1: Include hanging face nodes, but not edge nodes,\n"
      " h=2: Include all hanging nodes. (default)");
  // sc_options_add_string (sc_opts, 'v', "vtk", &opts->vtk, opts->vtk, "VTK basename");
  sc_options_add_string (sc_opts, 'm', "msh", &opts->mshpath, opts->mshpath, "MSH base filename");
  sc_options_add_string (sc_opts, 'c', "conn", &opts->conn, opts->conn, "Name of the connectivity");

  // sc_options_add_bool (sc_opts, 'h', "help", &help, help, "Print help message");
  // sc_options_add_bool (sc_opts, 0, "list-conn", &list_conn, list_conn, "Lists the available connectivity 'conn' options");

  exitcode = 0;
  do {
    first_argc = sc_options_parse(p4est_package_id, SC_LP_DEFAULT, sc_opts, argc, argv);
    if (first_argc < 0) {
      g->errmsg = "Failed to parse option format";
      sc_options_print_usage (p4est_package_id, SC_LP_ERROR, sc_opts, NULL);
      exitcode = 1;
      break;
    }

    if (first_argc != argc) {
      g->errmsg = "Too many arguments";
      sc_options_print_usage (p4est_package_id, SC_LP_ERROR, sc_opts, NULL);
      exitcode = 1;
      break;
    }

    if (help) {
      sc_options_print_usage (p4est_package_id, SC_LP_INFO, sc_opts, NULL);
      exitcode = 0;
      break;
    }

    if (list_conn) {
      char **conn_name = p4est_connectivity_names;
      for (;*conn_name; ++conn_name)
        printf("  %s\n", *conn_name);
      exitcode = 0;
      break;
    }

    g->opts = opts;
    // sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    if (create_context_from_opts(g, opts)) {
      exitcode = 1;
      break;
    }

    smesh = p4est_new_simplex_mesh(
        g->p4est,
        g->geom,
        g->ghost_layer, 
        opts->hangine_node_level);

#ifndef P4_TO_P8
    sprintf(filepath, "out2d/box_%s_%d-%d", opts->conn, opts->minlevel, opts->maxlevel);
#else
    sprintf(filepath, "out3d/box_%s_%d-%d", opts->conn, opts->minlevel, opts->maxlevel);
#endif

    p4est_vtk_write_file(g->p4est, g->geom, filepath);

    // sprintf(filepath, "out/sx_%s.msh", g->opts->conn);
    // write_gmsh_file(filepath, smesh);

#ifndef P4_TO_P8
    sprintf(filepath, "out2d/sx_%s_%d-%d_%d.vtk", opts->conn, opts->minlevel, opts->maxlevel, g->mpirank);
#else
    sprintf(filepath, "out3d/sx_%s_%d-%d_%d.vtk", opts->conn, opts->minlevel, opts->maxlevel, g->mpirank);
#endif
    p4est_simplex_mesh_write_vtk_file(filepath, smesh);

  } while (0);
  if (g->errmsg) {
    P4EST_GLOBAL_LERROR (g->errmsg);
  }

  if (smesh) {
    p4est_simplex_mesh_destroy(smesh);
  }

  /*** clean up and exit ***/
  destroy_context(g);

  sc_options_destroy (sc_opts);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return exitcode ? EXIT_FAILURE : EXIT_SUCCESS;
}


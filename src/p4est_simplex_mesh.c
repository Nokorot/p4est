// #include <p4est_to_p8est.h>

#ifndef P4_TO_P8
#include <p4est_nodes.h>
#include "p4est_simplex_mesh.h"
#else
#include <p8est_nodes.h>
#include "p8est_simplex_mesh.h"
#endif


// Local declarations

#ifndef P4_TO_P8
typedef struct {
    p4est_quadrant_t quad;
    p4est_locidx_t which_tree;
    int8_t face;

    p4est_locidx_t node_index;
}
p4est_quad_face_node_t;

unsigned int
p4est_quad_face_equal_fn(
    const p4est_quad_face_node_t *v1,
    const p4est_quad_face_node_t *v2,
    const void *used_data);

unsigned int
p4est_quad_face_hash_fn(const p4est_quad_face_node_t *v, const void *user_data);
#else 

typedef struct {
    p4est_quadrant_t quad;
    p4est_locidx_t which_tree;
    int8_t edge;

    p4est_locidx_t node_index;
}
p8est_quad_edge_node_t;

unsigned int
p8est_quad_edge_hash_fn(const p8est_quad_edge_node_t *v, const void *user_data);

unsigned int
p8est_quad_edge_equal_fn(
    const p8est_quad_edge_node_t *v1,
    const p8est_quad_edge_node_t *v2,
    const void *used_data);

#define p4est_record_hanging_data_t   p8est_record_hanging_data_t

#define p4est_record_hanging_faces    p8est_record_hanging_faces
#define p4est_record_per_quad         p8est_record_per_quad
#endif


typedef struct {
  p4est_t *p4est;
  p4est_geometry_t *geom;
  p4est_ghost_t *ghost;
  p4est_nodes_t *nodes;

  size_t num_indep;
  size_t num_fnodes;

  sc_hash_array_t *hash_record;
  sc_array_t *vertices;
} 
p4est_record_hanging_data_t;

sc_hash_array_t *
p4est_record_hanging_faces (
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost,
    p4est_nodes_t *nodes,
    sc_array_t *vertices);

void p4est_record_per_quad(
    p4est_record_hanging_data_t *d,
    p4est_topidx_t t,
    p4est_quadrant_t *quad,
    size_t elem,
    int is_ghost);


// Implementations

p4est_simplex_mesh_t *
p4est_new_simplex_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost,
    int hanging_node_level)
{
  p4est_simplex_mesh_t *mesh;
  p4est_nodes_t *nodes;
  sc_hash_array_t *hash_record;

  p4est_indep_t *cnode;

  size_t num_verticies;
  size_t num_elems, num_cnodes, num_fnodes;
  size_t off_elem_nodes, off_cnodes, off_fnodes;

#ifdef P4_TO_P8
  p8est_hang4_t *fnode;
  p8est_hang2_t *enode;

  size_t num_enodes;
  size_t off_enodes;
#else
  p4est_hang2_t *fnode;
#endif

#ifndef P4_TO_P8
  p4est_quad_face_node_t *found, qid;
#else
  p8est_quad_edge_node_t *found, qid;
#endif

  size_t zz, vi, fi, ei;
  p4est_locidx_t elem;
  p4est_topidx_t t, rt;
  p4est_tree_t *tree;
  sc_array_t *quadrants;
  p4est_quadrant_t *quad;
  
  int face;

  p4est_locidx_t cf, cc, c0, c1;

#ifdef P4_TO_P8
  int edge;
  int face_has_hanging_edge;

  p4est_locidx_t ce, c2, c3;

  p4est_locidx_t edge_noes[P8EST_EDGES];
#endif

  p4est_locidx_t *simplex;
  double abc[3], *vx;
  size_t position;

  // Generate p4est_nodes structure
  nodes = p4est_nodes_new(p4est, ghost);

  mesh = P4EST_ALLOC(p4est_simplex_mesh_t, 1);
  
  mesh->vertices = sc_array_new(3*sizeof(double));
  hash_record = p4est_record_hanging_faces(
                    p4est, geom, ghost, nodes, mesh->vertices);

    /// The vertex order is
    ///    ghost_nodes
    ///    cnodes
    ///    fnodes
    ///    enodes
    ///    elem_nodes
    ///    additional_face_nodes

  num_elems = nodes->num_local_quadrants;
  num_cnodes = nodes->indep_nodes.elem_count;
  num_fnodes = nodes->face_hangings.elem_count;

  off_cnodes = mesh->vertices->elem_count;
  off_fnodes = off_cnodes + num_cnodes;

#ifdef P4_TO_P8
  num_enodes = nodes->edge_hangings.elem_count;
  off_enodes = off_fnodes + num_fnodes;

  off_elem_nodes = off_enodes + num_enodes;

  num_verticies = num_elems + num_cnodes + num_fnodes + num_enodes;
#else
  off_elem_nodes = off_fnodes + num_fnodes;
  num_verticies = num_elems + num_cnodes + num_fnodes;
#endif

  sc_array_push_count(mesh->vertices, num_verticies);
  mesh->simplicies = sc_array_new((P4EST_DIM + 1) *sizeof(p4est_locidx_t));

  // Add all the vertices coming from independent nodes (indep_nodes)
  for (vi=0; vi < num_cnodes; ++vi) {
    cnode = sc_array_index(&nodes->indep_nodes, vi);
    vx = sc_array_index(mesh->vertices, off_cnodes + vi);

    abc[0] = (double) cnode->x / P4EST_ROOT_LEN;
    abc[1] = (double) cnode->y / P4EST_ROOT_LEN;
#ifdef P4_TO_P8
    abc[2] = (double) cnode->z / P4EST_ROOT_LEN;
#else
    abc[2] = 0;
#endif

    geom->X(geom, cnode->p.which_tree, abc, vx);
  }

  // Add all the vertices coming from face hangings nods (face_hangings)
  for (vi=0; vi < num_fnodes; ++vi) {
    fnode = sc_array_index(&nodes->face_hangings, vi);
    vx = sc_array_index(mesh->vertices, off_fnodes + vi);

    abc[0] = (double) fnode->x / P4EST_ROOT_LEN;
    abc[1] = (double) fnode->y / P4EST_ROOT_LEN;
#ifdef P4_TO_P8
    abc[2] = (double) fnode->z / P4EST_ROOT_LEN;
#else
    abc[2] = 0;
#endif

    geom->X(geom, fnode->p.which_tree, abc, vx);
  }

#ifdef P4_TO_P8
  // Add all the vertices coming from edge hangings nods (edge_hangings)
  for (vi=0; vi < num_enodes; ++vi) {
    enode = sc_array_index(&nodes->edge_hangings, vi);
    vx = sc_array_index(mesh->vertices, off_enodes + vi);

    abc[0] = (double) enode->x / P4EST_ROOT_LEN;
    abc[1] = (double) enode->y / P4EST_ROOT_LEN;
    abc[2] = (double) enode->z / P4EST_ROOT_LEN;

    geom->X(geom, enode->p.which_tree, abc, vx);
  }
#endif

  // We iterate through the local quads to make the simplices
  elem = 0;
  for (t=p4est->first_local_tree;
      t<=p4est->last_local_tree; ++t) {

    tree = p4est_tree_array_index(p4est->trees, t);
    qid.which_tree = t;

    quadrants = &tree->quadrants;
    for (zz = 0; zz < quadrants->elem_count; ++zz, ++elem) {
      quad = p4est_quadrant_array_index(quadrants, zz);
      qid.quad = *quad;

      // Add the quad centre node
      vx = sc_array_index(mesh->vertices, off_elem_nodes + elem);

      abc[0] = (double) ( quad->x + P4EST_QUADRANT_LEN(quad->level + 1) ) / P4EST_ROOT_LEN;
      abc[1] = (double) ( quad->y + P4EST_QUADRANT_LEN(quad->level + 1) ) / P4EST_ROOT_LEN;
#ifdef P4_TO_P8
      abc[2] = (double) ( quad->z + P4EST_QUADRANT_LEN(quad->level + 1) ) / P4EST_ROOT_LEN;
#else
      abc[2] = 0;
#endif
      geom->X(geom, t, abc, vx);
      cc = off_elem_nodes + elem;

#ifdef P4_TO_P8
      //  This should be in the record
      //  Could also be per face? This would avoid the first loop below
      //   However this way, there is only one lookup per quad in stead of per face

      for (ei=0; ei<P8EST_EDGES; ++ei) {
        qid.edge = ei;

        if (sc_hash_array_lookup(hash_record, &qid, &position)) {
          found = sc_array_index(&hash_record->a, position);

          if (found->node_index < 0)
            edge_noes[ei] = ~found->node_index;
          else 
            edge_noes[ei] = off_cnodes + found->node_index;
        }
        else
        {
          edge_noes[ei] = -1;
        }
      }
#endif

      for (fi = 0; fi < P4EST_FACES; ++fi) {
#ifndef P4_TO_P8
        qid.face = fi;

        c0 = off_cnodes +
          nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][0]];
        c1 = off_cnodes +
          nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][1]];

        // Check if the face is hanging
        if (sc_hash_array_lookup(hash_record, &qid, &position))
        {
          found = sc_array_index(&hash_record->a, position);
          if (found->node_index < 0)
            cf = ~found->node_index;
          else 
            cf = off_cnodes + found->node_index;

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = cf;

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = cf;
          simplex[2] = c1;
        }
        else
        {  // Split as usual
          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
        }
#else
        // Could skip this first loop if we know no edges are hanging
        face_has_hanging_edge = 0;
        for (ei = 0; ei < P4EST_HALF; ++ei) {
          edge = p8est_face_edges[fi][ei];
          if (0 <= edge_noes[edge])
          {
            face_has_hanging_edge = 1;
            break;
          }
        }

        if (face_has_hanging_edge)
        {
          cf = mesh->vertices->elem_count;
          vx = sc_array_push(mesh->vertices);

          // Compute face center vertex
          P4EST_ASSERT(0 <= fi && fi < 6);

          p4est_face_center_coords(quad, fi, abc);

          geom->X(geom, t, abc, vx);

          for (ei=0; ei<4; ++ei) {
            edge = p8est_face_edges[fi][ei];

            // Check these indices
            c0 = off_cnodes +
              nodes->local_nodes[P4EST_CHILDREN*elem + p8est_edge_corners[edge][0]];
            c1 = off_cnodes +
              nodes->local_nodes[P4EST_CHILDREN*elem + p8est_edge_corners[edge][1]];

            ce = edge_noes[edge];

            if (ce < 0)
            {
              simplex = sc_array_push(mesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = c1;
            }
            else
            {
              // edge center
              simplex = sc_array_push(mesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = ce;

              simplex = sc_array_push(mesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = ce;
              simplex[3] = c1;
            }
          }
        }
        else
        {
          c0 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][0]];
          c1 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][1]];
          c2 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][2]];
          c3 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][3]];

          int rot_corners = 0;
          if (quad->level == 0) {
            // The face split most be consistent on neighbouring faces.
            // This may be achieved by determining the diagonal, based on the global index 
            // of the vertices, since neighbouring faces has the same vertices. 
            // We do this by always splitting the diagonal with the vertex with the smallest 
            // global index.
            //
            // NOTE: We observe that for level = 0, all the vertices are independent in p4est_nodes_t
            // so we may assume they have a global index
    
            p4est_locidx_t first = c0;

            int rank = first < nodes->offset_owned_indeps 
                          ? nodes->nonlocal_ranks[first] 
                          : p4est->mpirank;
  

            // THIS IS SO HACKY!
            int rank1;
#define cmp_node(NODE) \
            rank1 = NODE < nodes->offset_owned_indeps  \
                          ? nodes->nonlocal_ranks[NODE] \
                          : p4est->mpirank; \
            if (rank1 < rank || (rank1 == rank && NODE < first)) { \
              first = NODE; \
              rank = rank1; \
            }

            cmp_node(c1);
            cmp_node(c2);
            cmp_node(c3);
#undef cmp_node

            rot_corners = (first == c1 || first == c2);
          } else {
            // We split the face, so that the parent face is split by an X:
            //    +---+---+
            //    | \ | / |
            //    +---+---+
            //    | / | \ |
            //    +---+---+
            // This way, the split is the same for all face-orientations (rotation, mirroring)
            // and thus the split will always be the same on any neighbouring face.
            //
            int dim = (fi >> 1);
            int c = p4est_quadrant_child_id(quad);
            int fchild =
                (dim == 0 ? c >> 1               : 0) 
              | (dim == 1 ? c & 1 | (c >> 1) & 2 : 0) 
              | (dim == 2 ? c & 3                : 0); 

            rot_corners = (fchild == 1 || fchild == 2);
          }

          // The face is (by default) split along the c0-c3 diagonal:
          //   c2--c3
          //   |  / |
          //   | /  |
          //   c0--c1
          //
          // In order to split along the opposite diagonal, we rotate the 
          // variables once anti-clockwise
          //
          //   c3--c1
          //   |  / |
          //   | /  |
          //   c2--c0
          //
          if (rot_corners) {
            p4est_locidx_t tmp;
            tmp = c1;
            c1 = c0;
            c0 = c2;
            c2 = c3;
            c3 = tmp;
          }

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
          simplex[3] = c3;

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c3;
          simplex[2] = c2;
          simplex[3] = c0;
        }
#endif
      }  // for (fi = 0..P4EST_FACES)
    } // for (zz = 0..|quadtrees|)
  } // for (t = first-..last-localtrees)

  sc_hash_array_destroy(hash_record);
  p4est_nodes_destroy(nodes);

  return mesh;
}

void
p4est_simplex_mesh_destroy(
    p4est_simplex_mesh_t *smesh)
{
  sc_array_destroy(smesh->vertices);
  sc_array_destroy(smesh->simplicies);
  P4EST_FREE(smesh);
}

void 
p4est_face_center_coords(
    p4est_quadrant_t *quad,
    int face,
    double abc[3])
{
  double frac = 1.0 / (double) P4EST_ROOT_LEN;
  double quad_len = P4EST_QUADRANT_LEN(quad->level);
  double half_len = P4EST_QUADRANT_LEN(quad->level + 1);

  abc[0] = frac * ( quad->x
           + ( !(face ^ 1) ? quad_len : 0)  // For face 1
           + (  (face >> 1) ? half_len : 0) // For all not x-orth faces
           );                             // (ie. fi != 0, 1)

  abc[1] = frac * ( quad->y
           + ( !(face ^ 3) ? quad_len : 0) // For face 3
           + ( !(face & 2) ? half_len : 0) // For all not y-orth faces
           );                            // (ie. fi != 2, 3)

#ifdef P4_TO_P8
  abc[2] = frac * ( quad->z
           + ( !(face ^ 5) ? quad_len : 0)  // For face 5
           + ( !(face >> 2) ? half_len : 0) // For all not z-orth faces
           );                             // (ie. fi != 4, 5)
#else
  abc[2] = 0;
#endif
}

#ifdef P4_TO_P8
void
p8est_edge_center_coords(
    p4est_quadrant_t *quad,
    int edge,
    double abc[3])
{
  double frac = 1.0 / (double) P4EST_ROOT_LEN;
  double len = P4EST_QUADRANT_LEN(quad->level);
  double half_len = 0.5 * len;
  
  int axis = (edge >> 2) + (edge >> 4);
  int num = (edge >> 2*axis) & 3;

  abc[0] = frac * ( quad->x
      + ( axis == 0 ? half_len : 0 )  // x-parallel
      + ( (axis != 0 && edge & 1) ? len : 0 )
      );
  abc[1] = frac * ( quad->y
      + ( axis == 1 ? half_len : 0 )  // y-parallel
      + ( (axis == 0 && edge & 1) || (axis == 2 && edge & 2) ? len : 0 )
      );
  abc[2] = frac * ( quad->z
      + ( axis == 2 ? half_len : 0 )  // z-parallel
      + ( (axis != 2 && edge & 2) ? len : 0 )
      );
}
#endif



#ifndef P4_TO_P8
unsigned int
p4est_quad_face_equal_fn(
    const p4est_quad_face_node_t *v1,
    const p4est_quad_face_node_t *v2,
    const void *used_data)
{
  return p4est_quadrant_equal_fn(&v1->quad, &v2->quad, NULL) \
    && v1->which_tree == v2->which_tree
    && v1->face == v2->face;
}

unsigned int
p4est_quad_face_hash_fn(const p4est_quad_face_node_t *v, const void *user_data)
{
  unsigned int qhash = p4est_quadrant_hash_fn(&v->quad, user_data);
  int b, c;
  b = v->which_tree;
  c = v->face;

  sc_hash_mix(qhash, b, c);
  sc_hash_final(qhash, b, c);

  return (unsigned) c;
}
#else 
unsigned int
p8est_quad_edge_hash_fn(const p8est_quad_edge_node_t *v, const void *user_data)
{
  unsigned int qhash = p4est_quadrant_hash_fn(&v->quad, user_data);
  int b, c;
  b = v->which_tree;
  c = v->edge;

  sc_hash_mix(qhash, b, c);
  sc_hash_final(qhash, b, c);

  return (unsigned) c;
}

unsigned int
p8est_quad_edge_equal_fn(
    const p8est_quad_edge_node_t *v1,
    const p8est_quad_edge_node_t *v2,
    const void *used_data)
{
  return p4est_quadrant_equal_fn(&v1->quad, &v2->quad, NULL) \
    && v1->which_tree == v2->which_tree
    && v1->edge == v2->edge;
}
#endif

void 
p4est_record_per_quad(
    p4est_record_hanging_data_t *d,
    p4est_topidx_t t,
    p4est_quadrant_t *quad,
    size_t elem,
    int is_ghost)
{
  int fi, face, nface;
  int c, fc, opp_corner0;

#ifdef P4_TO_P8
  int edge;
  int ec, opp_ec;

  p4est_locidx_t ce;
#endif

  double *vx, abc[3];

  p4est_topidx_t rt;
  p4est_quadrant_t parent, n;
  p4est_locidx_t c0;

#ifndef P4_TO_P8
  p4est_quad_face_node_t *added, qid;
#else
  p8est_quad_edge_node_t *added, qid;

  sc_array_t *quads, *treeids, *nedges;
#endif

  // if (quad->level == 0) { // This case need to be considered differently
  // TODO: Need to make sure faces are split the same from both sides
  //  Only relevant in 3D
  // }

  c = p4est_quadrant_child_id(quad);

  // Suffices to do this for child == 0, 7
  // if (c != 0 && c != P4EST_CHILDREN - 1)
  //   return;


  // Loop through the faces that are on the boundary of the parent,
  // and thus may be hanging
  for (fi = 0; fi < P4EST_DIM; ++fi) {
    face = p4est_corner_faces[c][fi];
    p4est_quadrant_parent(quad, &parent);

    if (is_ghost) {
        rt = p4est_quadrant_face_neighbor_extra(
                  &parent, t, face, &n, &nface, d->p4est->connectivity);
        if (rt < 0) // neighbour is outside conn
          continue;

        n.p.piggy3.which_tree = rt;

        ssize_t lnid; // p4est_quadrant_compare_piggy
        lnid = sc_array_bsearch (&(d->ghost->mirrors),
                                 &n,
                                 p4est_quadrant_compare_piggy);
        if (lnid < 0)  // neighbour \r n is not a local quad.
          continue;
    }
    else {
      fc = p4est_corner_face_corners[c][face];
      opp_corner0 = p4est_face_corners[face][(P4EST_HALF - 1) ^ fc];

      c0 = d->nodes->local_nodes[P4EST_CHILDREN * elem + opp_corner0];

      if (c0 < d->num_indep || c0 >= d->num_indep + d->num_fnodes)
        continue;

      rt = p4est_quadrant_face_neighbor_extra(
                &parent, t, face, &n, &nface, d->p4est->connectivity);

      if (rt < 0)
        continue;
    }

    // Hanging face
#ifndef P4_TO_P8
    qid.which_tree = rt;
    qid.quad = n;
    qid.face = nface % P4EST_FACES;

    added = sc_hash_array_insert_unique(d->hash_record, &qid, NULL);
    if (added) {

      *added = qid;
      if (is_ghost) {
        vx = sc_array_push(d->vertices);

        p4est_face_center_coords(&qid.quad, qid.face, abc);
        d->geom->X(d->geom, qid.which_tree, abc, vx);

        added->node_index = -d->vertices->elem_count;
      }
      else {
        added->node_index = c0;
      }
    }
#else
    for (size_t ei = 0; ei < 3; ++ei) {
      int edge = p8est_corner_edges[c][ei];
      int efe = p8est_edge_face_edges[edge][face];

      if (efe < 0)
        continue;

      int nedge;

      if (t == rt) {
        nedge = p8est_face_edges[nface][efe];
      } else {
        int nf = nface % P4EST_FACES, o = nface / P4EST_FACES;

        int nfedge =
            p8est_connectivity_face_neighbor_face_edge(efe, face, nf, o);

        nedge = p8est_face_edges[nf][nfedge];
      }

      qid.quad = n;
      qid.which_tree = rt;
      qid.edge = nedge;

      added = sc_hash_array_insert_unique(d->hash_record, &qid, NULL);
      if (added) {
        *added = qid;
  
        if (is_ghost) {
          vx = sc_array_push(d->vertices);
          ce = - d->vertices->elem_count;

          p8est_edge_center_coords(&qid.quad, qid.edge, abc);
          d->geom->X(d->geom, qid.which_tree, abc, vx);
        } else {
          ec = p8est_corner_edge_corners[c][edge];
          opp_ec = p8est_edge_corners[edge][1 ^ ec];

          ce = d->nodes->local_nodes[P4EST_CHILDREN * elem + opp_ec];
        }

        added->node_index = ce;
      }
    }
#endif
  }

#ifdef P4_TO_P8
  // Record hanging edges

  quads = sc_array_new(sizeof(p4est_quadrant_t));
  treeids = sc_array_new(sizeof(p4est_topidx_t));
  nedges = sc_array_new(sizeof(int));

  for (fi = 0; fi < P4EST_DIM; ++fi) {
    edge = p8est_corner_edges[c][fi];
    ec = p8est_edge_corners[edge][0];

    p4est_locidx_t c0 = d->nodes->local_nodes[P4EST_CHILDREN * elem + ec];

    if (c0 >= d->num_indep + d->num_fnodes) // Hanging edge
    {
      sc_array_reset(quads);
      sc_array_reset(treeids);
      sc_array_reset(nedges);

      p4est_quadrant_parent(quad, &parent);
      p8est_quadrant_edge_neighbor_extra(&parent, t, edge, quads, treeids,
                                         nedges, d->p4est->connectivity);

      for (int jj = 0; jj < quads->elem_count; ++jj) {
        qid.quad = *((p8est_quadrant_t *)sc_array_index(quads, jj));
        qid.which_tree = *((p4est_locidx_t *)sc_array_index(treeids, jj));
        qid.edge = *((int *)sc_array_index(nedges, jj));

        added = sc_hash_array_insert_unique(d->hash_record, &qid, NULL);
        if (added) {
          *added = qid;
          added->node_index = c0;
        }
      }
    }
  }

  sc_array_destroy(quads);
  sc_array_destroy(treeids);
  sc_array_destroy(nedges);
#endif
}

sc_hash_array_t *
p4est_record_hanging_faces (
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost,
    p4est_nodes_t *nodes,
    sc_array_t *vertices) 
{
  p4est_record_hanging_data_t d;

  p4est_topidx_t t, rt;
  p4est_tree_t *tree;
  p4est_quadrant_t *quad;

  sc_array_t *quadrants;
  size_t zz;
  
  d.p4est = p4est;
  d.geom = geom;
  d.ghost = ghost;
  d.nodes = nodes;
  d.vertices = vertices;
  d.num_indep = nodes->indep_nodes.elem_count;
  d.num_fnodes = nodes->face_hangings.elem_count;

#ifndef P4_TO_P8
  d.hash_record = sc_hash_array_new(sizeof(p4est_quad_face_node_t),
                                        (sc_hash_function_t) p4est_quad_face_hash_fn,
                                        (sc_equal_function_t) p4est_quad_face_equal_fn,
                                        NULL);
#else
  d.hash_record = sc_hash_array_new(sizeof(p8est_quad_edge_node_t),
                                        (sc_hash_function_t)  p8est_quad_edge_hash_fn,
                                        (sc_equal_function_t) p8est_quad_edge_equal_fn,
                                        NULL);
#endif

  p4est_locidx_t elem = 0;
  for (t=p4est->first_local_tree;
      t<=p4est->last_local_tree; ++t) {

    tree = p4est_tree_array_index(p4est->trees, t);
    quadrants = &tree->quadrants;

    for (zz = 0; zz < quadrants->elem_count; ++zz, ++elem) {
      quad = p4est_quadrant_array_index(quadrants, zz);
      p4est_record_per_quad(&d, t, quad, elem, 0);
    }
  }

  for (zz=0; zz < ghost->ghosts.elem_count; ++zz) {
    quad = sc_array_index(&(ghost->ghosts), zz);

    p4est_record_per_quad(&d, quad->p.which_tree, quad, 0, 1);
  }
  
  return d.hash_record;
}




/*** Export as gmsh/vtk file **/

void p4est_simplex_mesh_write_gmsh_file(const char *filename, p4est_simplex_mesh_t *mesh) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    // Write the header for Gmsh 2.2 format
    fprintf(file, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    // Write the nodes section
    size_t num_vertices = mesh->vertices->elem_count;
    fprintf(file, "$Nodes\n%zu\n", num_vertices);
    double *vertices = (double *) mesh->vertices->array;
    for (size_t i = 0; i < num_vertices; i++) {
        fprintf(file, "%zu %.16g %.16g %.16g\n", i + 1, vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
    }
    fprintf(file, "$EndNodes\n");

    // Write the elements section
    size_t num_simplicies = mesh->simplicies->elem_count;
    fprintf(file, "$Elements\n%zu\n", num_simplicies);
    p4est_locidx_t *simplicies = (p4est_locidx_t *) mesh->simplicies->array;
    for (size_t i = 0; i < num_simplicies; i++) {
#ifdef P4_TO_P8
        fprintf(file, "%zu 4 0 %u %u %u %u\n", i + 1, simplicies[4 * i] + 1, simplicies[4 * i + 1] + 1, simplicies[4 * i + 2] + 1, simplicies[4 * i + 3] + 1);
#else
        fprintf(file, "%zu 2 0 %u %u %u\n", i + 1, simplicies[3 * i] + 1, simplicies[3 * i + 1] + 1, simplicies[3 * i + 2] + 1);
#endif
    }
    fprintf(file, "$EndElements\n");

    fclose(file);
}

void p4est_simplex_mesh_write_vtk_file(const char *filename, p4est_simplex_mesh_t *mesh) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    // Write the header for VTK file
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Simplex Mesh\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    // Write the points section
    size_t num_vertices = mesh->vertices->elem_count;
    fprintf(file, "POINTS %zu double\n", num_vertices);
    double *vertices = (double *) mesh->vertices->array;
    for (size_t i = 0; i < num_vertices; i++) {
        fprintf(file, "%.16g %.16g %.16g\n", vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
    }

    // Write the cells section
    size_t num_simplicies = mesh->simplicies->elem_count;
#ifdef P4_TO_P8
    fprintf(file, "CELLS %zu %zu\n", num_simplicies, num_simplicies * 5);
#else
    fprintf(file, "CELLS %zu %zu\n", num_simplicies, num_simplicies * 4);
#endif
    p4est_locidx_t *simplicies = (p4est_locidx_t *) mesh->simplicies->array;
    for (size_t i = 0; i < num_simplicies; i++) {
#ifdef P4_TO_P8
        fprintf(file, "4 %u %u %u %u\n", simplicies[4 * i], simplicies[4 * i + 1], simplicies[4 * i + 2], simplicies[4 * i + 3]);
#else
        fprintf(file, "3 %u %u %u\n", simplicies[3 * i], simplicies[3 * i + 1], simplicies[3 * i + 2]);
#endif
    }

    // Write the cell types section
    fprintf(file, "CELL_TYPES %zu\n", num_simplicies);
    for (size_t i = 0; i < num_simplicies; i++) {
#ifdef P4_TO_P8
        fprintf(file, "10\n");  // VTK_TETRA type
#else
        fprintf(file, "5\n");   // VTK_TRIANGLE type
#endif
    }

    fclose(file);
}







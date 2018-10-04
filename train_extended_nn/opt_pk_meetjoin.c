
/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 */


/* ********************************************************************** */
/* opt_pk_meetjoin.c: Meet and join operations */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_test.h"
#include "opt_pk_meetjoin.h"
#include "opt_pk_project.h"
#include "opt_pk_cherni.h"
#include "rdtsc.h"
#include <time.h>
#include "python3.5m/Python.h"



unsigned short int step_size = 5;







typedef enum {
	STOER_MIN_CUT,
	UNARY_WEIGHT_DELETE,
	BINARY_WEIGHT_DELETE
}split_strategy;

typedef enum {
	NO_MERGE,
	MERGE_SMALL,
	MERGE_SMALL_LARGE
		
}merge_strategy;

#define max(a, b) (((a) > (b)) ? (a) : (b))



/* ====================================================================== */
/* II.2 Meet with (array of) linear constraint(s) */
/* ====================================================================== */

/* ---------------------------------------------------------------------- */
/* Factorized version */



bool opt_poly_meet_matrix(bool meet,
		      bool lazy,
		      elina_manager_t* man,
		      opt_pk_t* op,
		      opt_pk_t* oa, opt_matrix_t* oc)
{
  
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  man->result.flag_best = (oa->intdim==0);
  man->result.flag_exact = meet;
    size_t start = oa->C->nbrows;
    //assert(pa->satC);
    if(lazy){
	   
	    opt_poly_obtain_sorted_C(opk,oa);
	    if (op != oa){
	      op->C = opt_matrix_merge_sort(opk,oa->C,oc);
	      //op->nbeq = oa->nbeq;
	    }
	    else {
	        opt_matrix_merge_sort_with(opk,oa->C,oc);
	        if(oa->F){
			opt_matrix_free(oa->F);
			oa->F = NULL;
		}
		if(oa->satC){
			opt_satmat_free(oa->satC);
			oa->satC = NULL;
		}
		if(oa->satF){
			opt_satmat_free(oa->satF);
			oa->satF = NULL;
		}
		oa->nbline = 0;
		op->nbeq = oa->nbeq;
	    }
	oa->status = 0;
	return false;
     }
     else{
	size_t start = oa->C->nbrows;
    	assert(oa->satC);
		
	if (op != oa){
	      	op->C = opt_matrix_append(oa->C,oc);
	      	op->satC = opt_satmat_copy_resize_cols(oa->satC,
					opt_bitindex_size(op->C->nbrows));
	      	op->nbline = oa->nbline;
	      	op->nbeq = oa->nbeq;
	}
    	else {
      		opt_matrix_append_with(oa->C,oc);
      		opt_satmat_resize_cols(oa->satC,
			 opt_bitindex_size(oa->C->nbrows));
    	}
	opt_matrix_sort_rows(opk,op->F);
    	combine_satmat(opk,op,op->C->nbcolumns -2,start,true);
    	opt_cherni_add_and_minimize(opk,meet,op,start);
	if(opk->exn){
		
		return false;
	}
	if(!op->F){
		return true;
	}
	return false;
     }
    
  op->is_minimized = false;
  return false;
 
}




bool opt_poly_meet_elina_lincons_array(bool lazy,
				 elina_manager_t* man,
				 opt_pk_t* op, opt_pk_t* oa,
				 elina_lincons0_array_t* array)
{
  opt_matrix_t* mat;
  bool quasilinear;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  unsigned short int num_var = oa->intdim + oa->realdim;
  unsigned short int k; 
  quasilinear = elina_lincons0_array_is_quasilinear(array);
  /* quasilinearize if needed */
  if (!quasilinear){
    	elina_interval_t ** env = opt_generator_to_box(opk,oa->F);
    	quasilinearize_elina_lincons0_array(array,env,true,ELINA_SCALAR_MPQ);
	
	for(k=0; k < num_var; k++){
		elina_interval_free(env[k]);
	}
	free(env);
  }
  linearize_elina_lincons0_array(array,true,ELINA_SCALAR_MPQ);
  elina_lincons0_array_reduce_integer(array,op->intdim,ELINA_SCALAR_MPQ);
  bool exact = opt_matrix_set_elina_lincons0_array(opk,&mat,array,op->intdim,op->realdim,true);
  opt_matrix_sort_rows(opk,mat);
  size_t i;
  for(i=0; i < mat->nbrows; i++){
	if(!mat->p[i][0]){
		op->nbeq++;
	}
  }
  
  bool is_bottom = opt_poly_meet_matrix(true,lazy,man,op,oa,mat);
  
  opt_matrix_free(mat);
  if (opk->exn){
    	
    	man->result.flag_exact = man->result.flag_best = false;
  }
  else {
    	man->result.flag_best = man->result.flag_exact = exact ? true : false;
  }
  return is_bottom;
}

comp_list_t * linexpr0_to_comp_list(opt_pk_internal_t * opk, elina_linexpr0_t * expr){
	comp_list_t * cl = create_comp_list();
	size_t size = expr->size;
	size_t j;
	elina_linexpr_discr_t discr = expr->discr;
	if(discr==ELINA_LINEXPR_DENSE){
		elina_coeff_t* coeff = expr->p.coeff;
		for(j=0; j < size; j++){
			if(!elina_coeff_zero(&coeff[j])){
				insert_comp(cl,j + opk->dec);
			}
		}
	}
	else{
		elina_linterm_t* linterm = expr->p.linterm;
		for(j = 0; j < size; j++){
			elina_dim_t dim = linterm[j].dim;
			insert_comp(cl,dim + opk->dec);	
		}
	}
	return cl;
}

comp_list_t * lincons0_to_comp_list(opt_pk_internal_t * opk, elina_lincons0_t * cons){
	elina_linexpr0_t * expr = cons->linexpr0;
	return linexpr0_to_comp_list(opk,expr);
}

array_comp_list_t * lincons0_array_to_array_comp_list(opt_pk_internal_t *opk, elina_lincons0_array_t* array, unsigned short int n, char * is_trivial){
	size_t size = array->size;
	size_t i;
	array_comp_list_t * acl = create_array_comp_list();
	for(i = 0; i < size; i++){
		if(is_trivial[i]){
			continue;
		}
		elina_lincons0_t * cons = array->p + i;
		comp_list_t * cl = lincons0_to_comp_list(opk,cons);
		if(cl->size!=0){
			insert_comp_list_with_union(acl,cl,n);
		}			
	}
	return acl;
}

elina_linexpr0_t * copy_linexpr0_with_comp_list(opt_pk_internal_t *opk, elina_linexpr0_t * src, unsigned short int * ca, unsigned short int comp_size){
	elina_linexpr0_t * dst;
	elina_linexpr_discr_t discr = src->discr; 
	size_t size = src->size;
	unsigned short int j, k;
	if(discr==ELINA_LINEXPR_DENSE){
		dst = elina_linexpr0_alloc(discr,comp_size);
		elina_coeff_set(&dst->cst,&src->cst);
		elina_coeff_t* src_coeff = src->p.coeff;
		elina_coeff_t* dst_coeff = dst->p.coeff;
		for(k=0; k < comp_size; k++){
			unsigned short int num = ca[k] - opk->dec;
			if(!elina_coeff_zero(&src_coeff[num])){
				elina_coeff_set(&dst_coeff[k],&src_coeff[num]);
			}
		}
		
	}
	else{
		dst = elina_linexpr0_alloc(discr,size);
		elina_coeff_set(&dst->cst,&src->cst);
		elina_linterm_t* src_linterm = src->p.linterm;
		elina_linterm_t * dst_linterm = dst->p.linterm;
		
		for(j = 0; j < size; j++){
			elina_dim_t src_dim = src_linterm[j].dim + opk->dec;
			elina_coeff_t src_coeff = src_linterm[j].coeff;
			/************************
				Linexpr dimensions may not be sorted as observed in SeaHorn
			***********************/
			k = 0;
			while(ca[k] != src_dim){
				k++;
			}
			dst_linterm[j].dim = k;
			elina_coeff_set(&dst_linterm[j].coeff, &src_coeff);	
			k++;
		}
	}
	return dst;
}

void copy_lincons0_with_comp_list(opt_pk_internal_t *opk, elina_lincons0_t * dst, elina_lincons0_t * src, unsigned short int * ca, unsigned short int comp_size){
	dst->constyp = src->constyp;
        if(src->scalar!=NULL){
		if(dst->scalar==NULL){
			dst->scalar = elina_scalar_alloc();
		}
		elina_scalar_set(dst->scalar,src->scalar);
        }
	elina_linexpr0_t * src_expr = src->linexpr0;
	dst->linexpr0 = copy_linexpr0_with_comp_list(opk,src_expr,ca,comp_size);
}

bool is_linexpr_zero(elina_linexpr0_t * expr){
	elina_linexpr_discr_t discr = expr->discr;
	unsigned short int size = expr->size;
	unsigned short int j;	
	if(discr==ELINA_LINEXPR_DENSE){
		elina_coeff_t * coeff = expr->p.coeff;
		for(j=0; j < size; j++){
			if(!elina_coeff_zero(&coeff[j])){
				return false;
			}			
		}
	}
	else{
		elina_linterm_t * linterm = expr->p.linterm;
		for(j=0; j < size; j++){
			elina_coeff_t coeff = linterm[j].coeff;
			if(!elina_coeff_zero(&coeff)){
				return false;
			}
		}
	}
	return true;
}

int elina_coeff_sgn(elina_coeff_t * coeff){
    if(coeff->discr==ELINA_COEFF_SCALAR){
	return elina_scalar_sgn(coeff->val.scalar);
    }
    else{
	int si = elina_scalar_sgn(coeff->val.interval->inf);
	int ss = elina_scalar_sgn(coeff->val.interval->sup);
	if(!si){
		if(!ss){
			return 0;
		}
		else{
			return 2;
		}
	}
	else if(si > 0){
		return 1;
	}
	else{
		if(!ss){
			return -2;
		}
		else if(ss > 0){
			return -3;
		}
		else{
			return -1;
		}
	}
    }
	
}

opt_pk_array_t* opt_pk_meet_lincons_array_cons(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, elina_lincons0_array_t* array)
{
  //printf(".");
  //printf("meet lincons input\n");
  //elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,oa);
  //elina_lincons0_array_fprint(stdout,&arr2,NULL);
  //elina_lincons0_array_clear(&arr2);
  //elina_lincons0_array_fprint(stdout,array,NULL);
  //fflush(stdout);
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_MEET_LINCONS_ARRAY);
  size_t i;
  size_t size = array->size; 
  opt_pk_array_t * op;
  if(size==0){
	
	if(destructive){
		op = oa;
		return op;
	}
	else{
		op = opt_pk_copy(man,oa);
		return op;
	}
  }

  op = destructive ? oa :  opt_pk_array_alloc(NULL,NULL,oa->maxcols);
  array_comp_list_t * acla = oa->acl;
  
  /* if oa is bottom, return bottom */
  if(oa->is_bottom || !acla){
	 man->result.flag_best = man->result.flag_exact = true;
    	opt_poly_set_bottom(opk,op);
    	return op;
  }
  unsigned short int k;
  unsigned short int num_compa = acla->size;
  opt_pk_t ** poly_a = oa->poly;
  for(k=0; k < num_compa; k++){
      opt_pk_t * oak = poly_a[k];
	
      if(opk->funopt->algorithm>=0){
	 opt_poly_chernikova(man,oak,"meet lincons input");
      }
      else{
	 opt_poly_obtain_C(man,oak,"meet lincons input");
      }
      if(opk->exn){
	 opk->exn = ELINA_EXC_NONE;
	 if(!oak->C){
	    man->result.flag_best = man->result.flag_exact = false;
	    opt_poly_set_top(opk,op);
	    return op;
	 }
      }
      if(!oak->C && !oak->F){
	 man->result.flag_best = man->result.flag_exact = true;
   	 opt_poly_set_bottom(opk,op);
	 return op;
      }
  }
  bool is_bottom = false;
  char * is_trivial = (char *)calloc(size, sizeof(char));
  /*****************************************
	Handle Trivial constraints
  *****************************************/
  for(i=0; i < size; i++){
	elina_lincons0_t * cons = array->p + i;
	elina_linexpr0_t * expr = cons->linexpr0;
	elina_constyp_t constyp = cons->constyp;
	if(constyp==ELINA_CONS_DISEQ){
		is_trivial[i] = 1;
	}
	if(is_linexpr_zero(expr)){
		is_trivial[i] = 1;
		int sgn = elina_coeff_sgn(&expr->cst);
		if(constyp==ELINA_CONS_EQ){
			if(sgn){
				is_bottom = true;
				break;
			}
		}
		else if(constyp==ELINA_CONS_SUPEQ){
			if(sgn< 0){
				is_bottom = true;
				break;
			}
		}
		else if(constyp==ELINA_CONS_SUP){
			if(sgn != 1){
				is_bottom = true;
				break;
			}
		}
	}
  }

  if(is_bottom){
	man->result.flag_best = man->result.flag_exact = true;
    	opt_poly_set_bottom(opk,op);
    	return op;
  }


  array_comp_list_t * aclb = lincons0_array_to_array_comp_list(opk,array,oa->maxcols, is_trivial);
  
  /**********************************
	Compute Union of partitions
  **********************************/
  unsigned short int maxcols = oa->maxcols;
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
  unsigned short int num_comp = acl->size;
 
  unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
  unsigned short int * comp_size_map = (unsigned short int *)calloc(num_comp,sizeof(unsigned short int));
  comp_list_t * cl = acl->head;
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size; 
	ca_arr[k] = to_sorted_array(cl,maxcols);
	comp_size_map[k] = comp_size;
	cl = cl->next;
  }

  /**********************************
	Factor linear constraints according to union
  ***********************************/
  elina_lincons0_array_t * arr = (elina_lincons0_array_t *)malloc(num_comp*sizeof(elina_lincons0_array_t ));
  unsigned short int * rmapb = (unsigned short int *)calloc(size, sizeof(unsigned short int));
  size_t * nbmapb = (size_t *)calloc(num_comp,sizeof(size_t));
  for(i=0; i < size; i++){
	if(is_trivial[i]){
		continue;
	}
	elina_lincons0_t * cons = array->p + i;
	comp_list_t * clb = lincons0_to_comp_list(opk,cons);
	short int res = is_comp_list_included(acl,clb,maxcols);
	rmapb[i] = res;
        nbmapb[res]++;
	free_comp_list(clb);
  }
  
  for(k=0; k < num_comp; k++){
	arr[k] = elina_lincons0_array_make(nbmapb[k]);
  }
 
  size_t * counter = (size_t *)calloc(num_comp,sizeof(size_t));
  for(i =0; i < size; i++){
	if(is_trivial[i]){
		continue;
	}
	short int k1 = rmapb[i];
	unsigned short int * ca = ca_arr[k1];
        unsigned short int i1 = counter[k1];
	counter[k1]++;
	elina_lincons0_t * src = array->p + i;
	elina_lincons0_t * dst = (arr+k1)->p + i1;
	copy_lincons0_with_comp_list(opk,dst,src,ca, comp_size_map[k1]);
  }

  /*********************************
	Factor A according to union
  **********************************/
  
  unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));
  size_t * nbmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * nbeqmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  char * disjoint_map = (char *)calloc(num_compa,sizeof(char));
  size_t * nbgenmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * nblinemapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t * num_vertex = (size_t *)calloc(num_comp,sizeof(size_t));
  unsigned short int * array_map_a = create_array_map(acla,maxcols);
  comp_list_t * cla = acla->head;
  //unsigned short int * array_map_b = create_array_map(aclb,maxcols);
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	short int res = is_comp_list_included(acl,cla,maxcols);
	rmapa[k] = res;
	nbmapa[res] = nbmapa[res] + oak->C->nbrows;
	nbeqmapa[res] = nbeqmapa[res] + oak->nbeq;
	opt_poly_obtain_satF(oak);
	num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
	if(num_vertex_a[k]){
		if(!num_vertex[res]){
			num_vertex[res] = num_vertex_a[k];
		}
		else{
			num_vertex[res] = num_vertex[res] * num_vertex_a[k];
		}
	}
	nbgenmapa[res] = nbgenmapa[res] + oak->F->nbrows;
	nblinemapa[res] = nblinemapa[res] + oak->nbline;
	cla = cla->next;
  }
  
  opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  size_t * counterC = (size_t *)calloc(num_comp,sizeof(size_t));
  char * pos_con_map = (char *)calloc(num_comp,sizeof(char));
  size_t * counterF = (size_t *)calloc(num_comp,sizeof(size_t));
  cl = acl->head;
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size;
	poly[k] = opt_poly_alloc(comp_size,0);
	pos_con_map[k] = 1;
	cl = cl->next;
  }

  cl = acl->head;  	
  for(k=0; k < num_comp; k++){
	unsigned short int l,j;
	if(!nbmapb[k]){
		for(l=0; l < num_compa;l++){
			if(rmapa[l]==k){
				break;
			}
		}
		disjoint_map[l] = 1;
		if(destructive){
			poly[k]->C = poly_a[l]->C;
			poly[k]->nbeq = poly_a[l]->nbeq;
			poly[k]->F = poly_a[l]->F;
			poly[k]->satF = poly_a[l]->satF;
			poly[k]->satC = poly_a[l]->satC;
			poly[k]->nbline = poly_a[l]->nbline;
		}
		else{
			
			poly[k]->C = opt_matrix_copy(poly_a[l]->C);
			poly[k]->nbeq = poly_a[l]->nbeq;
			poly[k]->F = poly_a[l]->F ? opt_matrix_copy(poly_a[l]->F) : NULL; 
			poly[k]->satF = poly_a[l]->satF ? opt_satmat_copy(poly_a[l]->satF) : NULL; 
			poly[k]->satC = poly_a[l]->satC ? opt_satmat_copy(poly_a[l]->satC) : NULL; 
			poly[k]->nbline = poly_a[l]->nbline;
		}
	}
	else{
		unsigned short int comp_size = comp_size_map[k];
		poly[k]->C = opt_matrix_alloc(nbmapa[k]+1, comp_size+2,false);
		poly[k]->nbeq = nbeqmapa[k];		
		unsigned short int k1;
		unsigned short int nblines = 0;
		unsigned short int * ca = ca_arr[k];
		for(k1=0; k1 < comp_size; k1++){
			unsigned short int var = ca[k1];
			if(!array_map_a[var]){
				nblines++;
			}
		}
		poly[k]->F = opt_matrix_alloc(nbgenmapa[k]+2*num_vertex[k]+nblines+1, comp_size+2,false);
		num_vertex[k] = 0;
		nblines = 0;
		for(k1=0; k1 < comp_size; k1++){
			unsigned short int var = ca[k1];
			if(!array_map_a[var]){
				poly[k]->F->p[nblines][k1+2] = 1;
				nblines++;
			}
		}
		poly[k]->nbline = nblinemapa[k] + nblines;
		if(nblines==comp_size){
			poly[k]->F->p[nblines][0] = 1;
			poly[k]->F->p[nblines][1] = 1;
			nblines++;
		}
		counterF[k] = nblines;
		poly[k]->F->nbrows = nblines;
	}
	cl = cl->next;
  }
  
  free(array_map_a);
  // meet the constraints
  meet_cons_with_map(opk,oa,poly,rmapa, ca_arr, counterC, pos_con_map, disjoint_map);
  /*************************
		Add positivity constraint of A
  **************************/      
  for(k=0; k < num_comp; k++){
	size_t count = counterC[k];
	if(pos_con_map[k]&& nbmapb[k]){
		poly[k]->C->p[count][0] = 1;
		poly[k]->C->p[count][1] = 1;
		counterC[k]++;
		poly[k]->C->nbrows = counterC[k];
	}
  }	
	
  // cartesian product of vertices
  cartesian_product_vertices_with_map(oa, poly, rmapa, ca_arr, num_vertex_a, num_vertex, counterF, disjoint_map);

  // meet of rays
  meet_rays_with_map(oa, poly, rmapa, ca_arr, num_vertex_a, counterF, disjoint_map);
  char * exc_map = (char *)calloc(num_comp,sizeof(char));
  is_bottom = false;	
  for(k=0; k < num_comp; k++){
	if(nbmapb[k] && !is_bottom){
		poly[k]->satC = opt_satmat_alloc(poly[k]->F->nbrows,opt_bitindex_size(poly[k]->C->nbrows));
		is_bottom = opt_poly_meet_elina_lincons_array(opk->funopt->algorithm<0,
				      man,poly[k],poly[k],arr+k);
		if(opk->exn){
			opk->exn = ELINA_EXC_NONE;
			exc_map[k]=1;
		}
	  	
	}
	free(ca_arr[k]);
	elina_lincons0_array_clear(arr+k);
  }
  
  if(destructive){
	for(k=0; k < num_compa; k++){
		//unsigned short int ind = rmapa[k];
		if(!disjoint_map[k]){
			opt_poly_clear(poly_a[k]);
		}
		free(poly_a[k]);
	}
	free(poly_a);
  }
  
  if(!is_bottom){
	array_comp_list_t * tmpa = oa->acl;
	
	if(destructive){
		free_array_comp_list(tmpa);
	}
	unsigned short int k1=0;
	unsigned short int bound= num_comp;
	cl = acl->head;
	for(k=0; k < num_comp; k++){
	        opt_pk_t *oak = poly[k1];
	  	if(exc_map[k]){
	  		comp_list_t * tmp = cl;
	  		cl = cl->next;
	  		remove_comp_list(acl,tmp);
	  		unsigned short int k2;
	  		for(k2=k1; k2 < bound - 1; k2++){
	  			poly[k2] = poly[k2+1];
	  		}
	  		opt_poly_clear(oak);
	  		bound--;
	  	}
	  	else{
	  		k1++;
	  		cl=cl->next;
	  	}
	 }
	op->acl = acl;
	op->poly = poly;
	op->is_bottom = false;
  }
  else{
	for(k = 0; k < num_comp; k++){
		opt_pk_t * op_k = poly[k];
		if(op_k){
			opt_poly_clear(op_k);
			free(op_k);
		}
	}
	free(poly);
	free_array_comp_list(acl);
	op->is_bottom = true;
	
  }
  free(comp_size_map);
  free(rmapa);
  free(rmapb);
  free(nbmapa);
  free(nbmapb);
  free(nbeqmapa);
  free(counter);
  free(counterC);
  free(counterF);
  free(ca_arr);
  free(arr);
  free(pos_con_map);
  free(is_trivial);
  free(disjoint_map);
  free(num_vertex);
  free(num_vertex_a);
  free(nbgenmapa);
  free(nblinemapa);
  free(exc_map);
  free_array_comp_list(aclb);
  //for(k=0; k<num_comp; k++){
//	opt_matrix_fprint(stdout,poly[k]->F);
 // }
  //printf("meet lincons output\n");
  //elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,op);
  //elina_lincons0_array_fprint(stdout,&arr1,NULL);
  //elina_lincons0_array_clear(&arr1);
  	//fflush(stdout);
  return op;
}

opt_pk_array_t* opt_pk_meet_lincons_array(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, elina_lincons0_array_t* array){
	#if defined(TIMING)
		start_timing();
	#endif
	opt_pk_array_t * op = opt_pk_meet_lincons_array_cons(man,destructive,oa,array);
	#if defined(TIMING)
		record_timing(meet_lincons_time);
	#endif
	return op;
}
/************************************************
	Meet with tcons array
*************************************************/
opt_pk_array_t* opt_pk_meet_tcons_array(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, elina_tcons0_array_t* array)
{
  opt_pk_array_t *op;
  op= elina_generic_meet_intlinearize_tcons_array(man,destructive,oa,array,
						  ELINA_SCALAR_MPQ, ELINA_LINEXPR_LINEAR,
						  &opt_pk_meet_lincons_array);
  return op;
}



/* ********************************************************************** */
/* III. Join */
/* ********************************************************************** */

/* ====================================================================== */
/* III.1 Join of two or more polyhedra, functional and side-effect versions */
/* ====================================================================== */


bool are_vertex_disjoint(opt_pk_internal_t *opk, opt_numint_t * v1, opt_numint_t * v2, unsigned short int size){
	unsigned short int j;
	char * map1 = (char *)calloc(size - 2, sizeof(char));
	char * map2 = (char *)calloc(size - 2, sizeof(char));
	for(j = opk->dec; j < size; j++){
		map1[j-2] = (v1[j]> 0);
		map2[j-2] = (v2[j]> 0);
	}
	for(j=0; j < size-2; j++){
		if(map1[j]&&map2[j]){
			free(map1);
			free(map2);
			return false;
		}
	}
	free(map1);
	free(map2);
	return true;
}

bool is_vertex_included(opt_pk_internal_t *opk, opt_numint_t * v1, opt_numint_t * v2, unsigned short int size){
	unsigned short int j;
	char * map1 = (char *)calloc(size - 2, sizeof(char));
	char * map2 = (char *)calloc(size - 2, sizeof(char));
	for(j = opk->dec; j < size; j++){
		map1[j-2] = (v1[j]> 0);
		map2[j-2] = (v2[j]> 0);
	}
	for(j=0; j < size-2; j++){
		if(map1[j]&&map2[j]){
			if(v1[j+2]!=v2[j+2]){
				free(map1);
				free(map2);
				return false;
			}
		}
		if(map1[j]&&!map2[j]){
			free(map1);
			free(map2);
			return false;
		}
	}
	free(map1);
	free(map2);
	return true;
}


void combine_vertex(opt_pk_internal_t *opk, opt_numint_t *v1, opt_numint_t *v2, unsigned short int size){
	unsigned short int j;
	for(j=opk->dec; j < size; j++){
		v1[j] = v1[j] + v2[j];
	}
}


array_comp_list_t * compute_finest_partition(opt_matrix_t * C, comp_list_t *cl, unsigned short int dim, size_t * deleted_constraints){
	array_comp_list_t *res = create_array_comp_list();
	unsigned short int k=0;
	size_t i, nbrows = C->nbrows;
	unsigned short int nbcolumns = C->nbcolumns;
	unsigned short int *ca = to_sorted_array(cl,dim);
	unsigned short int j;
	for(i=0; i < nbrows; i++){
	    opt_numint_t *pi = C->p[i];
            if(deleted_constraints!=NULL && deleted_constraints[i]){
			continue;
	    }
	    comp_list_t * cl = create_comp_list();
	    for(j=2; j < nbcolumns; j++){
		if(pi[j]){
		   insert_comp(cl,ca[j-2]);
		}
	    }	
	    if(cl->size){
		  insert_comp_list_with_union(res,cl,dim);
	    }
        }
	free(ca);
	return res;
}


/***************************************
	Max Flow remove constraints
****************************************/

typedef struct node {
    int val;
    struct node *next;
} node_t;

void push(node_t *head, int val) {
    node_t *current = head;
    while (current->next != NULL) {
        current = current->next;
    }

    current->next = malloc(sizeof(node_t));
    current->next->val = val;
    current->next->next = NULL;
}

bool bfs(int **rGraph, int s, int t, int *parent, int size) {

    bool *visited = malloc(size * sizeof(bool));

    for (int i = 0; i < size; ++i)
        visited[i] = false;

    node_t *head = NULL;
    head = malloc(sizeof(node_t));
    if (head == NULL) {
        return 1;
    }

    head->val = s;
    head->next = NULL;

    visited[s] = true;
    parent[s] = -1;

    while (head != NULL) {
        int u = head->val;

        for (int v = 0; v < size; v++) {
            if (!visited[v] && rGraph[u][v] > 0) {
                push(head, v);
                parent[v] = u;
                visited[v] = true;
            }
        }
        head = head->next;
    }
    return (visited[t]);
}

void dfs(int **rGraph, int s, bool *visited, int size) {
    visited[s] = true;
    for (int i = 0; i < size; i++) {
        if (rGraph[s][i] > 0 && !visited[i])
            dfs(rGraph, i, visited, size);
    }
}

void partition_dfs(long long int **rGraph, int s, int *comp_num, int size, int last_comp) {
    comp_num[s] = last_comp;
    for (int i = 0; i < size; i++) {
        if (rGraph[s][i] > 0 && comp_num[i] == -1)
            partition_dfs(rGraph, i, comp_num, size, last_comp);
    }
}

int fordFulkerson(int **graph, int s, int t, int size, int maxFlow, int *x, int *y, bool *visited) {
    int u, v;

    int **rGraph = (int **) malloc(size * sizeof(int *));
    for (int i = 0; i < size; i++)
        rGraph[i] = (int *) malloc(size * sizeof(int));

    for (u = 0; u < size; u++)
        for (v = 0; v < size; v++)
            rGraph[u][v] = graph[u][v];


    int *parent = malloc(size * sizeof(int));
    for (int i = 0; i < size; ++i) {
        parent[i] = 0;
    }

    int max_flow = 0;  // There is no flow initially

    while (bfs(rGraph, s, t, parent, size)) {
        int path_flow = INT_MAX;
        for (v = t; v != s; v = parent[v]) {
//          printf("u=%d\n", u);
            u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }

        max_flow += path_flow;
    }


//  printf("\nmax-flow = %d , max=%d \n", max_flow, maxFlow);
    if (max_flow > maxFlow) {
//  printf("maxflow= %d", max_flow);

        for (int k = 0; k < size; ++k) {
            visited[k] = false;
        }

        dfs(rGraph, s, visited, size);
        x[0] = 1;
        y[0] = 1;
        for (int k = 1; k < size * size; ++k) {
            x[k] = -1;
            y[k] = -1;
        }
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                if (visited[i] && !visited[j] && graph[i][j] > 0) {
                    x[x[0]] = i;
                    x[0]++;
                    y[y[0]] = j;
                    y[0]++;
                }
    }
    for(int i =0; i < size; i++){
	free(rGraph[i]);
    }
    free(rGraph);
    free(parent);
    return max_flow;
}



/*******************************
	Stoer Wagner MinCut
********************************/
comp_list_t * stoer_wagner(long long int  **weights, unsigned short int N) 
{
    if(N==0){
	return NULL;
    }
    comp_list_t *cut = create_comp_list();
    comp_list_t * best_cut = NULL;
    char * used = (char *)calloc(N,sizeof(char));
    //for(int i=0; i < N; i++){
	//used[i] = N;
    //}
    //VI used(N), cut, best_cut;
    long long int best_weight = -1;
    char * added = (char *)calloc(N,sizeof(char));
    long long int * w = (long long int *)calloc(N,sizeof(long long int));
    for (int phase = N-1; phase >= 0; phase--) 
    {
        
	for(unsigned short int l =0; l < N; l++){
		w[l] = weights[0][l]; 
	} 
       
	for(unsigned short int l=0; l < N; l++){
		added[l] = used[l];
	}
        unsigned short int prev, last = 0;
        for (unsigned short int i = 0; i < phase; i++) 
        {
            prev = last;
            last = 65535;
            for (unsigned short int j = 1; j < N; j++)
            {
                if (!added[j] && (last == 65535 || w[j] > w[last])) {
			last = j;
		}
		
            }
		
            if (i == phase-1) 
            {
                for (unsigned short int j = 0; j < N; j++) weights[prev][j] += weights[last][j];
                for (unsigned short int j = 0; j < N; j++) weights[j][prev] = weights[prev][j];
                used[last] = 1;
		
		insert_comp(cut,last);
		
                //cut.push_back(last);
                if (best_weight == -1 || w[last] < best_weight) 
                {
		    if(best_cut!=NULL){	
		    	free_comp_list(best_cut);
		    }
		    best_cut = copy_comp_list(cut);
                    best_weight = w[last];
                }
            }
            else 
            {
                for (unsigned short int j = 0; j < N; j++)
                {
                    w[j] += weights[last][j];
                    added[last] = 1;
                }
            }
        }
    }
	free_comp_list(cut);
	free(used);
	free(added);
	free(w);
	unsigned short int best_cut_size = best_cut->size;
	if(best_cut_size < N/2){
		return best_cut;
	}
	else{
		char * map = create_map(best_cut,N);
		free_comp_list(best_cut);
		best_cut = create_comp_list();
		for(unsigned short int i=0; i < N; i++){
			if(!map[i]){
				insert_comp(best_cut,i);
			}
		}
		return best_cut;
	}
}

void **create_adj_matrix(opt_matrix_t * C, size_t * deleted_constraints, long long int ** graph, unsigned short int numvars, char * deleted_map) {
    opt_numint_t **mat = C->p;
    size_t nbrows = C->nbrows;
    unsigned short int nbcolumns = C->nbcolumns;
    for (unsigned short int i = 0; i < numvars; i++) {
	for(unsigned short int j=0; j < numvars; j++){
        	graph[i][j] = 0;
	}
    }

	
    for (size_t i = 0; i < nbrows; i++) {
	if(deleted_constraints[i]){
		continue;
	}
	unsigned short int j1 = 0;
        for (unsigned short int j = 2; j < nbcolumns; j++) {
	    
            if (mat[i][j] != 0) {
		unsigned short int k1 = 0;
                for (unsigned short int k = 2; k < j; k++) {
		    
                    if (mat[i][k] != 0) {
                        graph[j1][k1]++;
                        graph[k1][j1]++;
                    }
		    if(!deleted_map[k-2]){
			k1++;
	    	    }
                }
            }
	    if(!deleted_map[j-2]){
			j1++;
	    }
        }
    }
} 


array_comp_list_t * constraints_to_delete_stoer_wagner(opt_matrix_t *C, opt_matrix_t * newC, comp_list_t *cl, unsigned short int maxcols, unsigned short int threshold){
	
	size_t  nbrows = C->nbrows;
	unsigned short int nbvars = C->nbcolumns-2;
	size_t * deleted_constraints = (size_t *)calloc(nbrows,sizeof(size_t));
	long long int ** adj_mat = (long long int **)malloc(nbvars*sizeof(long long int*));
	for(unsigned short int j=0; j < nbvars; j++){
		adj_mat[j] = (long long int *)calloc(nbvars,sizeof(long long int));
	}
	comp_list_t * newcl = copy_comp_list(cl);
	opt_numint_t **p = C->p; 
	char * deleted_map = (char*)calloc(nbvars,sizeof(char));
	unsigned short int * ca_arr = to_sorted_array(cl,maxcols);	
	while(newcl->size > threshold){
		unsigned short int * ca = to_sorted_array(newcl,maxcols);
		unsigned short int * ind_map = map_index(ca,ca_arr,newcl->size);
        	create_adj_matrix(C, deleted_constraints,adj_mat,newcl->size,deleted_map);
		
		/*for(unsigned short int i=0; i < nbvars; i++){
			for(unsigned short int j = 0; j < nbvars; j++){
				printf("%u ",adj_mat[i][j]);
			}
			printf("\n");
		}*/
		comp_list_t * best_cut = stoer_wagner(adj_mat,newcl->size);
		unsigned short int best_cut_size = best_cut->size;
		
		
		for(size_t i =0; i < nbrows; i++){
			if(deleted_constraints[i]){
				continue;
			}
			opt_numint_t * pi = p[i];
			comp_t * c = best_cut->head;
			for(unsigned short int j=0; j < best_cut_size; j++){
				unsigned short int j1 = c->num;
				unsigned short int j2 = ind_map[j1];
				if(pi[j2+2]){
					deleted_constraints[i] = 1;
					//printf("deleted2: %d\n",i);
					break;
				}
				c = c->next;
			}
		}
		comp_t *c = best_cut->head;
		for(unsigned short int j=0; j < best_cut_size; j++){
			unsigned short int j1 = c->num;
			unsigned short int var = ca[j1];
			deleted_map[j1] = 1;
			remove_comp(newcl,var);	
			c = c->next;
		}
		
		//print_comp_list(newcl,maxcols);
		free_comp_list(best_cut);
		free(ca);
	}
	
	size_t i1=0;
	for(size_t i=0; i < nbrows; i++){
		if(!deleted_constraints[i]){
			opt_vector_copy(newC->p[i1],C->p[i],C->nbcolumns);
			i1++;	
		}
	}
	newC->nbrows = i1;
	array_comp_list_t * res = create_array_comp_list();
	insert_comp_list(res,newcl);
	
	
	for(unsigned short int j=0; j < nbvars; j++){
		free(adj_mat[j]);
	}
	free(adj_mat);
        free(deleted_constraints);
	//printf("coming here");
	//print_comp_list(cl,maxcols);
	//print_comp_list(newcl,maxcols);
	//opt_matrix_fprint(stdout,newC);
	//fflush(stdout);
	return res;
}

void merge_size_t(size_t *arr, size_t *arr2, int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    size_t L[n1], R[n2], L2[n1], R2[n2];

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) {
        L[i] = arr[l + i];
        L2[i] = arr2[l + i];
    }
    for (j = 0; j < n2; j++) {
	
        R[j] = arr[m + 1 + j];
        R2[j] = arr2[m + 1 + j];
    }

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
		
            arr[k] = L[i];
            arr2[k] = L2[i];
            i++;
        } else {
		
            arr[k] = R[j];
            arr2[k] = R2[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1) {
        arr[k] = L[i];
        arr2[k] = L2[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2) {
        arr[k] = R[j];
        arr2[k] = R2[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort_size_t(size_t *arr, size_t *arr2, int l, int r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort_size_t(arr, arr2, l, m);
        mergeSort_size_t(arr, arr2, m + 1, r);
	
        merge_size_t(arr, arr2, l, m, r);
	
    }
	
}


array_comp_list_t * constraints_to_delete(opt_matrix_t *C, opt_matrix_t * newC, comp_list_t *cl, unsigned short int maxcols, unsigned short int threshold){
	size_t  nbrows = C->nbrows;
	unsigned short int nbvars = C->nbcolumns-2;
	size_t * deleted_constraints = (size_t *)calloc(nbrows,sizeof(size_t));
	size_t * weights = (size_t *)calloc(nbvars,sizeof(size_t));
	size_t * sum_map = (size_t *)calloc(nbrows,sizeof(size_t));
	size_t * index_map = (size_t *)calloc(nbrows,sizeof(size_t));
	unsigned short int j;
	size_t i;
	for(i=0; i < nbrows; i++){
		for(j=0; j < nbvars; j++){
			if(C->p[i][j+2]!=0){
				weights[j]++;
			}
		}
	}
	
	
	for(i=0; i < nbrows; i++){
		size_t sum = 0;
		for(j=0; j < nbvars; j++){
			if(C->p[i][j+2]!=0){
				sum+=weights[j];
		        }
		}
		sum_map[i] = sum;
		index_map[i] = i;
	}
	array_comp_list_t * res = NULL;
	mergeSort_size_t(sum_map,index_map,0,(int)nbrows-1);
	for(i=0; i < nbrows; i++){
		size_t index = index_map[nbrows-1-i];
		
		deleted_constraints[index] = 1;
		array_comp_list_t *res = compute_finest_partition(C,cl,maxcols,deleted_constraints);
		
		bool flag = true;
		comp_list_t * cl = res->head;
		while(cl!=NULL){
			if(cl->size>threshold){
				flag = false;
			}
			cl = cl->next;
		}
		if(flag){
			size_t i1=0;
			for(i=0; i < nbrows; i++){
				if(!deleted_constraints[i]){
					
					opt_vector_copy(newC->p[i1],C->p[i],C->nbcolumns);
					i1++;	
				}
			}
			newC->nbrows = i1;
			free(weights);
			free(sum_map);
			free(index_map);
			free(deleted_constraints);
			return res;	
		}
	}
	free(weights);
	free(sum_map);
	free(index_map);
        free(deleted_constraints);
	size_t i1=0;
	for(i=0; i < nbrows; i++){
		if(!deleted_constraints[i]){
			opt_vector_copy(newC->p[i1],C->p[i],C->nbcolumns);
			i1++;	
		}
	}
	newC->nbrows = i1;
	return res;
}


array_comp_list_t * constraints_to_delete_binary(opt_matrix_t *C, opt_matrix_t * newC, comp_list_t *cl, unsigned short int maxcols, unsigned short int threshold){
	size_t  nbrows = C->nbrows;
	unsigned short int nbvars = C->nbcolumns-2;
	size_t * deleted_constraints = (size_t *)calloc(nbrows,sizeof(size_t));
	size_t ** graph = (size_t **)calloc(nbvars,sizeof(size_t*));
	size_t * sum_map = (size_t *)calloc(nbrows,sizeof(size_t));
	size_t * index_map = (size_t *)calloc(nbrows,sizeof(size_t));
	unsigned short int j,k;
	size_t i;
	for(j=0; j < nbvars; j++){
		graph[j] = (size_t*)calloc(nbvars,sizeof(size_t));
	}

	for(i=0; i < nbrows; i++){
		for(j=0; j < nbvars; j++){
			if(C->p[i][j+2]!=0){
				for(k=j+1; k<nbvars; k++){
					if(C->p[i][k+2]!=0){
						graph[j][k]++;
						graph[k][j]++;
					}
				}
			}
		}
	}
	
	
	for(i=0; i < nbrows; i++){
		size_t sum = 0;
		for(j=0; j < nbvars; j++){
			if(C->p[i][j+2]!=0){
				for(k=j+1; k < nbvars; k++){
					if(C->p[i][k+2]!=0){
						sum+=graph[j][k];
					}
				}
		        }
		}
		sum_map[i] = sum;
		index_map[i] = i;
	}
	array_comp_list_t * res = NULL;
	mergeSort_size_t(sum_map,index_map,0,(int)nbrows-1);
	for(i=0; i < nbrows; i++){
		size_t index = index_map[nbrows-1-i];
		
		deleted_constraints[index] = 1;
		array_comp_list_t *res = compute_finest_partition(C,cl,maxcols,deleted_constraints);
		
		bool flag = true;
		comp_list_t * cl = res->head;
		while(cl!=NULL){
			if(cl->size>threshold){
				flag = false;
			}
			cl = cl->next;
		}
		if(flag){
			size_t i1=0;
			for(i=0; i < nbrows; i++){
				if(!deleted_constraints[i]){
					
					opt_vector_copy(newC->p[i1],C->p[i],C->nbcolumns);
					i1++;	
				}
			}
			newC->nbrows = i1;
			for(j=0; j < nbvars; j++){
				free(graph[j]);
			}
			free(graph);
			free(sum_map);
			free(index_map);
			free(deleted_constraints);
			return res;	
		}
	}
	for(j=0; j < nbvars; j++){
		free(graph[j]);
	}
	free(graph);
	free(sum_map);
	free(index_map);
        free(deleted_constraints);
	size_t i1=0;
	for(i=0; i < nbrows; i++){
		if(!deleted_constraints[i]){
			opt_vector_copy(newC->p[i1],C->p[i],C->nbcolumns);
			i1++;	
		}
	}
	newC->nbrows = i1;
	return res;
}


int *find_deleted_constraints(int size, int consSize, long long int **cMat) {

    int **graph = (int **) malloc(size * sizeof(int *));
    int **seenMap = (int **) malloc(size * sizeof(int *));
    for (int i = 0; i < size; i++) {
        graph[i] = (int *) malloc(size * sizeof(int));
        seenMap[i] = (int *) malloc(size* sizeof(int));
    }

    int *s = malloc(size * size * sizeof(int));
    int *t = malloc(size * size * sizeof(int));

    for (int k = 0; k < size * size; ++k) {
        s[k] = -1;
        t[k] = -1;
    }

    s[0] = 1;
    t[0] = 1;
   /* size_t * weight = (size_t *)calloc(size-2,sizeof(size_t));
    for(int i = 0 ; i < consSize; i++){
	for(int j=2; j < size; j++){
		if(cMat[i][j]!=0){
			weight[j-2]++;
		}
	}
    }*/
     
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            graph[i][j] = 0;
	    seenMap[i][j] = 0;
        }
    }
	printf("size: %d\n",size);
	fflush(stdout);
    for (int i = 0; i < consSize; i++) {
        for (int j = 2; j < size; j++) {
            if (cMat[i][j] != 0) {
                for (int k = j+1; k < size; k++) {
			
                    if (cMat[i][k] != 0) {
                        graph[j - 2][k - 2]++;//consSize-weight[j-2] - weight[k-2];
                        graph[k - 2][j - 2]++;//consSize - weight[j-2] - weight[k-2];
                        if (!seenMap[j - 2][k - 2] && !seenMap[k - 2][j - 2]) {
                            seenMap[j - 2][k - 2] = 1;
                            seenMap[k - 2][j - 2] = 1;
                            s[s[0]] = j - 2;
                            t[t[0]] = k - 2;
                            s[0]++;
                            t[0]++;
                        }

                    }
                }
            }
        }
    }

    printf("graph\n");
	fflush(stdout);
    for(int i=0; i <size; i++){
	for(int j=0; j < size; j++){
		printf("%d ", graph[i][j]);
	}
	printf("\n");
    }
	
	 printf("s\n");
	fflush(stdout);
    for(int i=0; i <size; i++){
	for(int j=0; j < size; j++){
		printf("%d ", s[size*i+j]);
	}
	printf("\n");
    }
	
	 printf("t\n");
	fflush(stdout);
    for(int i=0; i <size; i++){
	for(int j=0; j < size; j++){
		printf("%d ", t[i*size+j]);
	}
	printf("\n");
    }
	fflush(stdout);
    int maxFlow = 0;
    int *x = malloc(size * size * sizeof(int));
    int *y = malloc(size * size * sizeof(int));

    for (int k = 0; k < size * size; ++k) {
        x[k] = -1;
        y[k] = -1;
    }
    x[0] = 1;
    y[0] = 1;
    bool *visited = malloc(size * sizeof(bool));

    for (int i = 1; i < size * size && s[i] != -1; i++) {
        maxFlow = max(fordFulkerson(graph, s[i], t[i], size, maxFlow, x, y, visited), maxFlow);
    }
	 /*printf("x\n");
	fflush(stdout);
    for(int i=0; i <size; i++){
	for(int j=0; j < size; j++){
		printf("%d ", x[size*i+j]);
	}
	printf("\n");
    }
	
	 printf("y\n");
	fflush(stdout);
    for(int i=0; i <size; i++){
	for(int j=0; j < size; j++){
		printf("%d ", y[i*size+j]);
	}
	printf("\n");
    }
	fflush(stdout);*/
    int *deletedConstraints = malloc(consSize * sizeof(int));
    for (int i = 0; i < consSize; ++i) {
        deletedConstraints[i] = 0;
    }
//    printf("printing visited vector: %p\n", visited);
//    fflush(stdout);
//    for (int i = 0; i < size; ++i) {
//        printf("%d\t", visited[i]);
//    }
//    printf("\n");
//    fflush(stdout);

    for (int k = 0; k < size * size && x[k] != -1; k++) {
        int xIndex = x[k] + 2;
        int yIndex = y[k] + 2;
	printf("x: %d y: %d\n",x[k],y[k]);
	fflush(stdout);
        for (int i = 0; i < consSize; i++) {
            if (cMat[i][xIndex] != 0 && cMat[i][yIndex] != 0) {
                deletedConstraints[i] = 1;
            }
        }
    }
  for(int i=0; i < size; i++){
	free(graph[i]);
	free(seenMap[i]);
  }
  free(graph);
  free(seenMap);
  free(s);
  free(t);
  free(x);
  free(y);
  free(visited);
//    for (int i = 0; i < consSize; i++) {
//        printf("\n%d ", deletedConstraints[i]);
//        fflush(stdout);
//    }
    return deletedConstraints;
}


/*void find_number_of_partitions(opt_matrix_t *C) {

    long long int **graph = create_adj_matrix(C->p, C->nbcolumns - 2, C->nbrows);

    printf("\nprinting g\n");
    fflush(stdout);
    for (int i = 0; i < C->nbcolumns - 2; ++i) {
        for (int j = 0; j < C->nbcolumns - 2; ++j) {
            printf("%lld \t ", graph[i][j]);
        }
        printf("\n");
        fflush(stdout);
    }
    printf("done\n");
    fflush(stdout);

    int *partition_map = malloc((C->nbcolumns - 2) * sizeof(int));
    for (int i = 0; i < C->nbcolumns - 2; ++i) {
        partition_map[i] = -1;
    }
    int last_partition = 0;

    for (int j = 0; j < C->nbcolumns - 2; ++j) {
        if (partition_map[j] == -1) {
            partition_dfs(graph, j, partition_map, C->nbcolumns - 2, last_partition);
            last_partition++;
        }

    }
    printf("\nprinting partition map\n");
    fflush(stdout);
    for (int l = 0; l < C->nbcolumns - 2; ++l) {
        printf("%d\t", partition_map[l]);
        fflush(stdout);
    }
    printf("done\n");
    fflush(stdout);

}*/


comp_list_t *
eliminate_constraint_rows(opt_pk_t *poly, opt_matrix_t *oldC, size_t *deletedConstraints, unsigned short int *sortedVars,
                          opt_pk_internal_t *opk) {
    size_t newNbRows = 0;
    printf("original\n");
    opt_matrix_fprint(stdout,oldC);
    fflush(stdout);
    for (size_t i = 0; i < poly->C->nbrows; i++) {
        if (deletedConstraints[i] == 1) newNbRows++;
    }
    newNbRows = poly->C->nbrows - newNbRows + 1;

    long long int **newP = (long long int **) malloc(newNbRows * sizeof(long long int *));
    for (size_t i = 0; i < newNbRows; i++) {
        newP[i] = (long long int *) malloc(poly->C->nbcolumns * sizeof(long long int));
    }
    for (size_t k = 0; k < newNbRows; ++k) {
        for (int i = 0; i < poly->C->nbcolumns; ++i) {
            newP[k][i] = 0;
        }
    }
	
    
    size_t counter = 0;
    for (size_t i = 0; i < poly->C->nbrows; i++) {
        if (deletedConstraints[i] == 0) {
            for (size_t j = 0; j < poly->C->nbcolumns; j++) {
                newP[counter][j] = poly->C->p[i][j];
            }
            counter++;
        }
        else{
		opt_vector_print(poly->C->p[i],oldC->nbcolumns);
	}
    }
	
    newP[counter][0] = 1;
    newP[counter][1] = 1;

    comp_list_t *newCl = create_comp_list();

    //find new number of variables
    size_t new_comp_num = 0;

    for (size_t i = 2; i < poly->C->nbcolumns; ++i) {
        for (size_t j = 0; j < newNbRows; ++j) {
            if (newP[j][i] != 0) {
                insert_comp(newCl, sortedVars[i - 2]);
                new_comp_num++;
                break;
            }
        }
    }


    unsigned short int *map_idx = (unsigned short int *) calloc(new_comp_num, sizeof(unsigned short int));
    size_t last_idx = 0;

    for (size_t i = 2; i < poly->C->nbcolumns; ++i) {
        for (size_t j = 0; j < newNbRows; ++j) {
            if (newP[j][i] != 0) {
                map_idx[last_idx] = (unsigned short) (i - 2);
                last_idx++;
                break;
            }
        }
    }


    oldC->nbrows = (size_t) newNbRows;
    oldC->_maxrows = (size_t) newNbRows;
    oldC->p = newP;


    if (new_comp_num == 0) {
        oldC = poly->C;
        return NULL;
    }
    split_matrix(opk, oldC, opt_matrix_copy(oldC), map_idx, new_comp_num);
    oldC->nbcolumns = (unsigned short) (new_comp_num + 2);
//    find_number_of_partitions(C);

    return newCl;
}


// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(unsigned short int *arr, unsigned short int *arr2, unsigned short int l, unsigned short int m, unsigned short int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    int L[n1], R[n2], L2[n1], R2[n2];

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) {
        L[i] = arr[l + i];
        L2[i] = arr2[l + i];
    }
    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1 + j];
        R2[j] = arr2[m + 1 + j];
    }

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            arr2[k] = L2[i];
            i++;
        } else {
            arr[k] = R[j];
            arr2[k] = R2[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1) {
        arr[k] = L[i];
        arr2[k] = L2[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2) {
        arr[k] = R[j];
        arr2[k] = R2[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(unsigned short int *arr, unsigned short int *arr2, unsigned short int l, unsigned short int r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr, arr2, l, m);
        mergeSort(arr, arr2, m + 1, r);

        merge(arr, arr2, l, m, r);
    }
}


comp_list_t *find_comp_by_index(comp_list_t *head, int num) {
    int i = 0;
    comp_list_t *res = head;
    while (i < num) {
        i++;
        res = res->next;
    }
    return res;
}

void create_joined_partition(elina_manager_t *man, array_comp_list_t *res, opt_pk_t **poly, opt_pk_t **poly1,
                             opt_pk_t **poly2,
                             array_comp_list_t *acl,
                             comp_list_t *clp, size_t nbhF, size_t nbhC, size_t nbheq, size_t nbhline,
                             size_t num_vertex1ha,
                             size_t num_vertex2hb, size_t nbhlinea, size_t *num_vertex1, size_t *num_vertex2,
                             unsigned short int maxcols, size_t *nblinemapa, unsigned short int **ca_arr,
                             char *pos_con_map, opt_pk_internal_t *opk, int index, unsigned short int mapFlag, unsigned short int num_comp_res) {

	
    unsigned short int comp_size = clp->size;
    poly[index] = opt_poly_alloc(comp_size, 0);
    poly[index]->F = opt_matrix_alloc(nbhF + 2 * (num_vertex1ha + num_vertex2hb), comp_size + 2, false);
    poly[index]->C = opt_matrix_alloc(nbhC, comp_size + 2, false);
    poly[index]->nbeq = nbheq;
    poly[index]->nbline = nbhline;
    unsigned short int *ca = to_sorted_array(clp, maxcols);
    size_t count = 0;
    //combine all overlapping vertices of A into one
    count = cartesian_product_vertices_one_comp(poly1, acl, ca_arr, num_vertex1, poly[index], count, ca,
                                                    pos_con_map, mapFlag);
	
    //combine all overlapping rays of A into one
    meet_rays_one_comp(poly1, acl, ca_arr, nblinemapa, poly[index], num_vertex1, ca, count, pos_con_map, mapFlag);
    size_t begin = poly[index]->F->nbrows - nbhlinea;
    //combine all overlapping constraints of A into one
    bool is_pos_con = meet_cons_one_comp(opk, poly1, acl, ca_arr, poly[index], ca, pos_con_map, mapFlag);
	
    if (is_pos_con) {
        size_t count = poly[index]->C->nbrows;
        poly[index]->C->p[count][0] = 1;
        poly[index]->C->p[count][1] = 1;
        count++;
        poly[index]->C->nbrows = count;
    }
	
    //combine all overlapping vertices of B into one
    count = cartesian_product_vertices_one_comp(poly2, acl, ca_arr, num_vertex2, poly[index],
                                                    poly[index]->F->nbrows, ca, pos_con_map, mapFlag);
    //combine all overlapping rays of B into one
    meet_rays_one_comp(poly2, acl, ca_arr, NULL, poly[index], num_vertex2, ca, count, pos_con_map, mapFlag);
    free(ca);
    opt_matrix_t *F = poly[index]->F;
    opt_matrix_t *C = poly[index]->C;
    //remove_common_gen(opk,F,begin);
    opt_matrix_sort_rows(opk, C);
    /************************
        Combine satmat of A
    *************************/
	
    poly[index]->satF = opt_satmat_alloc(C->nbrows, opt_bitindex_size(F->nbrows));
    combine_satmat(opk, poly[index], comp_size, begin, false);
    size_t num = F->nbrows - begin;
    opt_matrix_sort_rows_from(opk, F, begin, num);
    opt_poly_dual(poly[index]);
    opt_cherni_add_and_minimize(opk, false, poly[index], begin);
    opt_poly_dual(poly[index]);
	
    if(opk->exn){
	opt_pk_t *tmp = poly[index];
	unsigned short int k1;
	for(k1=index; k1 < num_comp_res-1; k1++){
		poly[k1] = poly[k1+1];
	}
	opt_poly_clear(tmp);
	//opk->exn= ELINA_EXC_NONE;
    }
    else{
    	insert_comp_list(res, copy_comp_list(clp));
    }
	
}

/*
void merge_blocks_large_small(elina_manager_t *man, opt_pk_internal_t *opk, opt_pk_t **poly, opt_pk_t **poly1, opt_pk_t **poly2,
			      size_t * num_vertex1, size_t * num_vertex2, size_t * nblinemapa, 
			      array_comp_list_t *res, array_comp_list_t *acl, unsigned short int **ca_arr, char * pos_con_map,
			      unsigned short int * join_index, unsigned short int maxcols, unsigned short int counter,
			      unsigned short int k2, unsigned short int nbcomp){
		unsigned short int left = 0, right = counter - 1;
        	while (left <= right) {
            		char *clp_map = (char *) calloc(maxcols, sizeof(char));
            		comp_list_t *clp = create_comp_list();
            		size_t num_vertex1a = 0, num_vertex2b = 0;
            		size_t nbF = 0, nbC = 0, nbeq = 0, nbline = 0, nblinea = 0;
            		if (left == right) {
                		unsigned short int kp = join_index[left];
                		union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
				if (!num_vertex1a) {
				    num_vertex1a = num_vertex1[kp];
				} else {
				    num_vertex1a *= num_vertex1[kp];
				}
				if (!num_vertex2b) {
				    num_vertex2b = num_vertex2[kp];
				} else {
				    num_vertex2b *= num_vertex2[kp];
				}
				nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
				nbC = nbC + poly1[kp]->C->nbrows;
				nbeq = nbeq + poly1[kp]->nbeq;
				nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
				nblinea += nblinemapa[kp];
                		create_joined_partition(man, res, poly, poly1, poly2, acl, clp, nbF, nbC, nbeq, nbline, num_vertex1a,
                                        num_vertex2b,
                                        nblinea, num_vertex1, num_vertex2, maxcols, nblinemapa, ca_arr, pos_con_map,
                                        opk, k2, pos_con_map[kp], nbcomp);
				if(opk->exn==ELINA_EXC_OVERFLOW){
					opk->exn = ELINA_EXC_NONE;
					nbcomp--;
				}
			
                		break;
            		}
            		unsigned short int sum_p = join_size[right];
            		unsigned short int tmp_sum = join_size[right];
            		unsigned short int kp = join_index[right];
            		union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
		        if (!num_vertex1a) {
		        	num_vertex1a = num_vertex1[kp];
		    	} else {
		        	num_vertex1a *= num_vertex1[kp];
		    	}
		    	if (!num_vertex2b) {
		        	num_vertex2b = num_vertex2[kp];
		    	} else {
		        	num_vertex2b *= num_vertex2[kp];
		    	}
            		nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
			
            		nbC = nbC + poly1[kp]->C->nbrows;
            		nbeq = nbeq + poly1[kp]->nbeq;
            		nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
            		nblinea += nblinemapa[kp];

            		while (left < right) {
                		tmp_sum += join_size[left];
                		if (tmp_sum > threshold) break;

                		sum_p += join_size[left];
                		kp = join_index[left];
                		union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
				if (!num_vertex1a) {
				    num_vertex1a = num_vertex1[kp];
				} else {
				    num_vertex1a *= num_vertex1[kp];
				}
				if (!num_vertex2b) {
				    num_vertex2b = num_vertex2[kp];
				} else {
				    num_vertex2b *= num_vertex2[kp];
				}
				nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
				
				nbC = nbC + poly1[kp]->C->nbrows;
				nbeq = nbeq + poly1[kp]->nbeq;
				nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
				nblinea += nblinemapa[kp];
				left++;
            		}
			
            		create_joined_partition(man, res, poly, poly1, poly2, acl, clp, nbF, nbC, nbeq, nbline, num_vertex1a,
                                    num_vertex2b,
                                    nblinea, num_vertex1, num_vertex2, maxcols, nblinemapa, ca_arr, pos_con_map, opk,
                                    k2, pos_con_map[kp], nbcomp);
			if(opk->exn==ELINA_EXC_OVERFLOW){
				opk->exn = ELINA_EXC_NONE;
				nbcomp--;
			}
			
            		k2--;
           		right--;
			free(clp_map);
}*/



/*****************************
		Join based on generator representation
*******************************/
opt_pk_array_t * opt_poly_join_gen(elina_manager_t *man, opt_pk_array_t *oa, opt_pk_array_t *ob, bool destructive){
	//printf("JOIN INPUT\n");
	//elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa);
        //elina_lincons0_array_fprint(stdout,&arr1,NULL);
	//elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,ob);
        //elina_lincons0_array_fprint(stdout,&arr2,NULL);
	//elina_lincons0_array_clear(&arr1);
	//elina_lincons0_array_clear(&arr2);
	//print_array_comp_list(oa->acl,oa->maxcols);
	//print_array_comp_list(ob->acl,ob->maxcols);
	//fflush(stdout);
	srand ( time(NULL) );
	unsigned short int num_thresholds = 0;
  	for(unsigned short int threshold = step_size; threshold < oa->maxcols;threshold=threshold+step_size){
		num_thresholds++;
  	}
	opt_pk_array_t * input1 = opt_pk_copy(man,oa);
	join_count++;
	opt_pk_internal_t *opk = opt_pk_init_from_manager(man, ELINA_FUNID_JOIN);
	size_t i;
	unsigned short int j,k;
	unsigned short int maxcols = oa->maxcols;
	array_comp_list_t * acla = oa->acl;
	array_comp_list_t * aclb = ob->acl;
	opt_pk_array_t *op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,maxcols);
	if(oa->is_bottom || !acla){
		if(destructive){
			opt_poly_array_clear(opk,op);
			free(op);
		}
		else{
			free(op);
		}
		op = opt_pk_copy(man,ob);
		
		return op;
	}
	
	if(ob->is_bottom || !aclb){
		if(destructive){
			op = oa;
		}
		else{
			free(op);
			op = opt_pk_copy(man,oa);
		}
		
		return op;
	}
	unsigned short int * map_a = create_array_map(acla,maxcols);
        unsigned short int * map_b = create_array_map(aclb,maxcols);
	elina_dim_t * tdim_a = (elina_dim_t *)calloc(maxcols,sizeof(elina_dim_t));
	size_t size_a = 0;
	comp_list_t *cla = acla->head;
	while(cla!=NULL){
		unsigned short int *ca = to_sorted_array(cla,maxcols);
		unsigned short int comp_size = cla->size;
		for(k=0; k < comp_size; k++){
			unsigned short int num = ca[k];
			if(!map_b[num]){
				tdim_a[size_a] = num-opk->dec;
				size_a++;
			}
		}
		free(ca);
		cla = cla->next;
	}
	
	opt_pk_array_t * tmp_a = NULL;
	if(size_a){
		if(destructive){
			tmp_a = opt_pk_forget_array(man,true,oa,tdim_a,size_a,false);
		}
		else{
			tmp_a = opt_pk_forget_array(man,false,oa,tdim_a,size_a,false);
		}
	}
	else{
		 tmp_a = oa;
	}	
	
	
	elina_dim_t * tdim_b = (elina_dim_t *)calloc(maxcols,sizeof(elina_dim_t));
	size_t size_b = 0;
	comp_list_t *clb = aclb->head;
	while(clb!=NULL){
		unsigned short int *ca = to_sorted_array(clb,maxcols);
		unsigned short int comp_size = clb->size;
		for(k=0; k < comp_size; k++){
			unsigned short int num = ca[k];
			if(!map_a[num]){
				tdim_b[size_b] = num-opk->dec;
				size_b++;
			}
		}
		free(ca);
		clb = clb->next;
	}

	opt_pk_array_t * tmp_b = NULL;
	if(size_b){
		tmp_b = opt_pk_forget_array(man,false,ob,tdim_b,size_b, false);
	}
	else{
		tmp_b = ob;
	}
	free(map_b);
	free(map_a);
	free(tdim_a);	
	free(tdim_b);

	/***********************
			Minimize A and compute components that could be connected by sigmas
	************************/
	opt_pk_t **poly_a = tmp_a->poly;
	acla = tmp_a->acl;
	unsigned short int num_compa = acla->size;
	cla = acla->head;
	size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
	for(k=0; k < num_compa; k++){
		opt_pk_t * oak = poly_a[k];
		
		if(opk->funopt->algorithm>=0){
			opt_poly_chernikova(man,oak,"cons to gen");
		}
		else{	
			opt_poly_obtain_F(man,oak,"cons to gen");
		}
		if(opk->exn){
			opt_poly_set_top(opk,op);
			
			return op;
		}
		if(!oak->F){
			if(destructive){
				opt_poly_array_clear(opk,op);	
  				free(op);
			}
			else{
				free(op);
			}
			op = opt_pk_copy(man,ob);
			return op;
		}
		/*if(destructive){
			if(oak->C){
				opt_matrix_free(oak->C);
				oak->C = NULL;
			}
			if(oak->satC){
				opt_satmat_free(oak->satC);
				oak->satC = NULL;
			}
			if(oak->satF){
				opt_satmat_free(oak->satF);
				oak->satF = NULL;
			}
		}*/
		//printf("A \n");
		//opt_matrix_fprint(stdout,oak->F);
		//opt_matrix_fprint(stdout,oak->C);
		//opt_satmat_fprint(stdout,oak->satF);
		//fflush(stdout);
		opt_poly_obtain_satF(oak);
		num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
		cla = cla->next;
	}	
	
	
	/***********************
			Minimize B and compute components that could be connected by sigmas
	************************/
	opt_pk_t **poly_b = tmp_b->poly;
	aclb = tmp_b->acl;
	unsigned short int num_compb = aclb->size;
	clb = aclb->head;
	size_t * num_vertex_b = (size_t *)calloc(num_compb,sizeof(size_t));
	for(k=0; k < num_compb; k++){
		opt_pk_t * obk = poly_b[k];
		if(opk->funopt->algorithm >=0){
			opt_poly_chernikova(man,obk,"cons to gen");
		}
		else{
			opt_poly_obtain_F(man,obk,"cons to gen");
		}
		if(opk->exn){
			opt_poly_set_top(opk,op);
			
			return op;
		}
		if(!obk->F){
			if(destructive){
				return oa;
			}
			else{
				free(op);
				op = opt_pk_copy(man,oa);
				return op;
			}
		}
		opt_poly_obtain_satF(obk);
		num_vertex_b[k] = opt_generator_rearrange(obk->F,obk->satF);
		clb = clb->next;
	}
	
	/**************************
			Compute union of independent components
	***************************/
	array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
	
	unsigned short int num_comp = acl->size;
	char ** var_map_a = (char **)malloc(num_comp*sizeof(char*));
	char ** var_map_b = (char **)malloc(num_comp*sizeof(char*));
	opt_pk_t ** poly1 = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
 	opt_pk_t ** poly2 = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
	size_t * num_vertex1 = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * num_vertex2 = (size_t *)calloc(num_comp,sizeof(size_t));	
	size_t * nbmapCa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nblinemapa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbeqmapa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbmapCb = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nblinemapb = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbeqmapb = (size_t *)calloc(num_comp,sizeof(size_t));		
	unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));


	comp_list_t * cl = acl->head;
        for(k=0; k < num_comp; k++){
		var_map_a[k] = (char *)calloc(cl->size, sizeof(char)); 
		var_map_b[k] = (char *)calloc(cl->size, sizeof(char)); 
		ca_arr[k] = to_sorted_array(cl,maxcols);
		cl = cl->next;
	}
	/**************************
		Partition A according to union 
	***************************/
	unsigned short int * rmapa = (unsigned short int *)calloc(num_compa,sizeof(unsigned short int));
	size_t * nbmapFa = (size_t *)calloc(num_comp,sizeof(size_t));
	
	cla = acla->head;
	for(k=0; k < num_compa; k++){
		opt_pk_t *oak = poly_a[k];
		unsigned short int inda = is_comp_list_included(acl,cla,maxcols);
		rmapa[k] = inda;	
		if(num_vertex_a[k]){
			if(!num_vertex1[inda]){
				num_vertex1[inda] = num_vertex_a[k];
			}
			else{
				num_vertex1[inda] = num_vertex1[inda] * num_vertex_a[k];
			}
		}
		
		nbmapFa[inda] = nbmapFa[inda] + oak->F->nbrows; 
		nbmapCa[inda] = nbmapCa[inda] + oak->C->nbrows; 
		nblinemapa[inda] = nblinemapa[inda] + oak->nbline;
		nbeqmapa[inda] = nbeqmapa[inda] + oak->nbeq;
		unsigned short int *ca = to_sorted_array(cla,maxcols);
		unsigned short int * ind_map = map_index(ca,ca_arr[inda],cla->size);
		unsigned short int k1;
		for(k1=0; k1 < cla->size; k1++){
			unsigned short int num = ind_map[k1];
			var_map_a[inda][num]++;
		}
		free(ca);
		free(ind_map);
		cla = cla->next;
	}
 
	/**************************
		Partition B according to union
	***************************/
	unsigned short int * rmapb = (unsigned short int *)calloc(num_compb,sizeof(unsigned short int));
	
	size_t * nbmapFb = (size_t *)calloc(num_comp,sizeof(size_t));
	clb = aclb->head;
	for(k=0; k < num_compb; k++){
		opt_pk_t * obk = poly_b[k];
		unsigned short int indb = is_comp_list_included(acl,clb,maxcols);
		rmapb[k] = indb;
		
		
		if(num_vertex_b[k]){
			if(!num_vertex2[indb]){
				num_vertex2[indb] = num_vertex_b[k];
			}
			else{
				num_vertex2[indb] = num_vertex2[indb] * num_vertex_b[k];
			}
		}
		nbmapCb[indb] = nbmapCb[indb] + obk->C->nbrows;
		nbeqmapb[indb] = nbeqmapb[indb] + obk->nbeq;
		nbmapFb[indb] = nbmapFb[indb] + obk->F->nbrows;
		nblinemapb[indb] = nblinemapb[indb] + obk->nbline;
		//printf("B %d\n",num_vertex_b[k]);
		//opt_matrix_fprint(stdout,obk->F);
		//opt_matrix_fprint(stdout,obk->C);
		//opt_satmat_fprint(stdout,obk->satF);
		//fflush(stdout);
		unsigned short int *ca = to_sorted_array(clb,maxcols);
		unsigned short int * ind_map = map_index(ca,ca_arr[indb],clb->size);
		unsigned short int k1;
		for(k1=0; k1 < clb->size; k1++){
			unsigned short int num = ind_map[k1];
			var_map_b[indb][num]++;
		}
		free(ca);
		free(ind_map);
		clb = clb->next;
	}
	
	cl = acl->head;
	size_t * counterFa = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterFb = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterCa = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterCb = (size_t *)calloc(num_comp, sizeof(size_t));
	//char ** vertex_map = (char **)malloc(num_comp*sizeof(char *));
	for(k = 0; k < num_comp; k++){
		unsigned short int comp_size = cl->size;
		size_t gen_size, nblines=0, nbeq=0;
		unsigned short int k1;
		
		// A
		poly1[k] = opt_poly_alloc(comp_size,0);
		if(!nbmapFa[k]){
		    gen_size = 1;
		}
		else{
			nblines = 0;			
			for(k1=0; k1 < cl->size; k1++){
				if(!var_map_a[k][k1]){
					nblines++;
					nbeq++;
				}
			}
			gen_size = nbmapFa[k] + 2*num_vertex1[k] + nblines;
		}
		poly1[k]->F = opt_matrix_alloc(gen_size,comp_size+2,false);
		poly1[k]->C = opt_matrix_alloc(nbmapCa[k]+1+nbeq,comp_size+2,false);
		poly1[k]->nbline = nblinemapa[k]+nblines;
		poly1[k]->nbeq = nbeqmapa[k]+nbeq;
		
		// B
		poly2[k] = opt_poly_alloc(comp_size,0);
		if(!nbmapFb[k]){
		    gen_size =  1;
		}
		else{
		    nblines = 0;
		    for(k1=0; k1 < cl->size; k1++){
			if(!var_map_b[k][k1]){
				nblines++;
			}
		    }
		    gen_size = nbmapFb[k] + 2*num_vertex2[k] + nblines;
		}
		
		poly2[k]->F = opt_matrix_alloc(gen_size,comp_size+2,false);
		poly2[k]->C = opt_matrix_alloc(nbmapCb[k]+1,comp_size+2,false);
		num_vertex1[k] = 0;
		num_vertex2[k] = 0;
		poly2[k]->nbline = nblinemapb[k]+nblines;
		poly2[k]->nbeq = nbeqmapb[k];
		cl = cl->next;
		
	}
	
	
	/*************************
		Cartesian Product of Vertices from  A
	************************/
	
	cartesian_product_vertices(tmp_a,poly1,rmapa,ca_arr,num_vertex_a,num_vertex1,counterFa);
	
	/************************
		Now Consider rays of A
	************************/
	meet_rays(tmp_a,poly1,rmapa,ca_arr,num_vertex_a,counterFa);
	
	//for(k=0; k < num_comp; k++){
	//	size_t count = counterF[k];
	//	start_counter[k] = count;
	//}

	/***********************
		Meet constraints of A
	************************/
	char * pos_con_map = (char *)calloc(num_comp, sizeof(char)); 
	for(k=0; k < num_comp; k++){
		pos_con_map[k] = 1;
	}
	meet_cons(opk,tmp_a,poly1,rmapa,ca_arr,counterCa,pos_con_map);
		
	/*************************
		Add positivity constraint of A
	**************************/
        cl = acl->head;
	for(k=0; k < num_comp; k++){
		size_t countC = counterCa[k];
		size_t countF = counterFa[k];
		nblinemapa[k] = 0;
		if(pos_con_map[k]){
			poly1[k]->C->p[countC][0] = 1;
			poly1[k]->C->p[countC][1] = 1;
			countC++;
		}
		if(nbmapFa[k]){
			unsigned short int k1;
			for(k1=0; k1 < cl->size; k1++){
				if(!var_map_a[k][k1]){
					poly1[k]->F->p[countF][k1+2]=1;
					poly1[k]->C->p[countC][k1+2]=1;
					nblinemapa[k]++;
					countC++;
					countF++;
				}
			}
		}
		counterCa[k] = countC;
		counterFa[k] = countF;
		poly1[k]->C->nbrows = counterCa[k];
		poly1[k]->F->nbrows = counterFa[k];
		cl = cl->next;
	}	
	
	
	
	/*************************
		Cartesian Product of Vertices from B
	************************/
	cartesian_product_vertices(tmp_b,poly2,rmapb,ca_arr,num_vertex_b,num_vertex2,counterFb);

	/************************
		 Consider rays of B
	************************/
	meet_rays(tmp_b,poly2,rmapb,ca_arr,num_vertex_b,counterFb);	
	
	
	
	/***********************
		Meet constraints of B
	************************/ 
	for(k=0; k < num_comp; k++){
		pos_con_map[k] = 1;
	}
	meet_cons(opk,tmp_b,poly2,rmapb,ca_arr,counterCb,pos_con_map);
	
	/*************************
		Add positivity constraint of B
	**************************/
        cl = acl->head;
	for(k=0; k < num_comp; k++){
		size_t countC = counterCb[k];
		size_t countF = counterFb[k];
		if(pos_con_map[k]){
			poly2[k]->C->p[countC][0] = 1;
			poly2[k]->C->p[countC][1] = 1;
			countC++;
		}
		if(nbmapFb[k]){
			unsigned short int k1;
			for(k1=0; k1 < cl->size; k1++){
				if(!var_map_b[k][k1]){
					poly2[k]->F->p[countF][k1+2]=1;
					countF++;
				}
			}
		}
		counterCb[k] = countC;
		counterFb[k] = countF;
		poly2[k]->C->nbrows = counterCb[k];
		poly2[k]->F->nbrows = counterFb[k];
		cl = cl->next;
	}	
	

	unsigned short int num_comp_res = 0;
	bool flag = false;
	cl = acl->head;
	unsigned short int counter = 0;
	unsigned short int * join_index = (unsigned short int *)calloc(num_comp,sizeof(unsigned short int));
	unsigned short int * join_size = (unsigned short int *)calloc(num_comp,sizeof(unsigned short int));
	
    unsigned short int min_block=maxcols, max_block=0, num_blocks=0, sum_blocks =0;
  	
	size_t max_gen = 0, min_gen = 18446744073709551615UL,  sum_gen = 0;
	double avg_gen = 0, avg_blocks = 0;
	
	elina_interval_t ** interval_a= opt_pk_to_box(man,oa);
	elina_interval_t ** interval_b = opt_pk_to_box(man,ob);
	for(k=0; k < num_comp; k++){
		
		pos_con_map[k] = 0;
		opt_matrix_t * Ca = poly1[k]->C;
		opt_matrix_t * Fa = poly1[k]->F;
		
		opt_matrix_t * Cb = poly2[k]->C;
		opt_matrix_t * Fb = poly2[k]->F;
		
		size_t num_gen = Fa->nbrows + Fb->nbrows;
		if(!nbmapFa[k] || !nbmapFb[k]){
			pos_con_map[k] = 3;
		}
		else if(opt_poly_leq(opk,Ca,Fb)&&opt_poly_leq(opk,Cb,Fa)){
			num_comp_res++;
			// take B
			pos_con_map[k] = 2;
		}
		else{
			num_blocks++;
				
			sum_blocks = sum_blocks + cl->size;
			sum_gen = sum_gen + num_gen;
			if(cl->size < min_block){
				min_block = cl->size;
			}
			if(cl->size > max_block){
				max_block = cl->size;
			}
			if(num_gen < min_gen){
				min_gen = num_gen;
			}
			if(num_gen > max_gen){
				max_gen = num_gen;
			}
		    		
			
		}
		cl = cl->next;
	}
	
	unsigned short int threshold = maxcols;

	if(max_block <= 0){
		min_gen = 0;
		avg_gen = 0;
		avg_blocks = 0;
		max_gen = 0;
		min_block = 0;
		max_block = 0;
	}
  	int ss = rand()%3;
  	int ms = rand()%3;

	int num_inf = 0;
	int num_semi_inf = 0;
	int num_bounded = 0;
	int num_exact = 0;
	
	cl = acl->head;
	char * reward_map = (char *)calloc(maxcols-2,sizeof(char));
	for(k=0; k < num_comp; k++){
		if(pos_con_map[k]==2|| pos_con_map[k]==3){
			cl = cl->next;
			continue;
		}
		else if(cl->size>threshold){
			
			num_comp_res++;
			pos_con_map[k] = 4;
			comp_t * c = cl->head;
			while(c!=NULL){
				unsigned short int num = c->num-opk->dec;
				reward_map[num] = 1;
				elina_interval_t * itv_a = interval_a[num];
				elina_interval_t * itv_b = interval_b[num];

				if(elina_interval_is_top(itv_a)){
					num_inf++;
				} else if(!elina_scalar_infty(itv_a->sup) && !elina_scalar_infty(itv_a->inf)){
					if(elina_scalar_equal(itv_a->sup,itv_a->inf)){
						num_exact++;
					}else{
						num_bounded++;
					}
				} else {
					num_semi_inf++;
				}
				
				if(elina_interval_is_top(itv_b)){
					num_inf++;
				} else if(!elina_scalar_infty(itv_b->sup) && !elina_scalar_infty(itv_b->inf)){
					if(elina_scalar_equal(itv_b->sup,itv_b->inf)){
						num_exact++;
					}else{
						num_bounded++;
					}
				} else {
					num_semi_inf++;
				}
				c = c->next;
			}
		}
		else{
			
			join_size[counter] = cl->size;
			join_index[counter] = k;
			counter++;
			comp_t * c = cl->head;
			while(c!=NULL){
				unsigned short int num = c->num-opk->dec;
				reward_map[num] = 1;
				elina_interval_t * itv_a = interval_a[num];
				elina_interval_t * itv_b = interval_b[num];
			
				if(elina_interval_is_top(itv_a)){
					num_inf++;
				} else if(!elina_scalar_infty(itv_a->sup) && !elina_scalar_infty(itv_a->inf)){
					if(elina_scalar_equal(itv_a->sup,itv_a->inf)){
						num_exact++;
					}else{
						num_bounded++;
					}
				} else {
					num_semi_inf++;
				}
				
				if(elina_interval_is_top(itv_b)){
					num_inf++;
				} else if(!elina_scalar_infty(itv_b->sup) && !elina_scalar_infty(itv_b->inf)){
					if(elina_scalar_equal(itv_b->sup,itv_b->inf)){
						num_exact++;
					}else{
						num_bounded++;
					}
				} else {
					num_semi_inf++;
				}
				c = c->next;
			}
			//flag = true;
		}
		cl =cl->next;
	}

	for(unsigned short int k3 = 0; k3 < oa->maxcols - 2; k3++){
		elina_interval_free(interval_a[k3]);
		elina_interval_free(interval_b[k3]);
	}
	free(interval_a);
	free(interval_b);

	// PyEval_AcquireLock();

	// if (!PyEval_ThreadsInitialized())
	// {
	// 	PyEval_InitThreads();
	// }
	// PyInterpreterState* Plock = PyInterpreterState_New();
	// PyEval_RestoreThread(Plock);

	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();

	PyObject* pName = PyUnicode_DecodeFSDefault("call_nn_double");
	PyObject* pModule = PyImport_Import(pName);
	
	bool max_block_gt_0 = false;
	if(max_block > 0){
		max_block_gt_0 = true;

		avg_blocks = ((double)sum_blocks)/num_blocks;
		avg_gen = ((double)sum_gen)/num_blocks;

		unsigned short int num_features = 11;
		int state[num_features];
		state[0] = num_blocks;
		state[1] = min_block; 
		state[2] = max_block;
		state[3] = avg_blocks;
		state[4] = min_gen; 
		state[5] = max_gen;
		state[6] = avg_gen;
		state[7] = num_inf;
		state[8] = num_semi_inf;
		state[9] = num_bounded;
		state[10] = num_exact;

		int action = 0;


		PyObject* pFunc = PyObject_GetAttrString(pModule, "predict");


		PyObject* pArgs = PyTuple_New(num_features);
		for(int i = 0; i < num_features; i++){
			PyObject* pValue = PyLong_FromLong(state[i]);
			PyTuple_SetItem(pArgs, i, pValue);
		}
		
		PyObject* res = PyObject_CallObject(pFunc, pArgs);
		action = (int) PyLong_AsLong(res);
		// printf("the nn guesses : %i\n",action);

		// printf("%i\n",action);
		short int tmp;
		threshold = (action/9+1)*step_size;
		tmp = action%9;
		ss = tmp/3;
		ms = tmp%3;

		
	}

	clock_t started = clock();

	if(counter>=1){
		 mergeSort(join_size, join_index, 0, counter - 1);
		 unsigned short int free_map = 5;
        	 unsigned short int left = 0, right = counter - 1;
        	 unsigned short int map_idx = 0;
		if(ms==NO_MERGE){
			while(left<=right){
				map_idx = join_index[left];
				pos_con_map[map_idx] = free_map;
				free_map++;
				num_comp_res++;
				left++;
			}
		}
		else{
			 while (left <= right) {
		    		if (left == right) {
				    map_idx = join_index[right];
		        	    pos_con_map[map_idx] = free_map;
		        	    free_map++;
		        	    num_comp_res++;
		        	    break;
		    		}
				unsigned short int  tmp_sum = 0;
				if(ms==MERGE_SMALL_LARGE){
					map_idx = join_index[right];
					tmp_sum = join_size[right];
					pos_con_map[map_idx] = free_map;
		    			num_comp_res++;
					right--;
				}
				else{
					map_idx = join_index[left];
					pos_con_map[map_idx] = free_map;
					num_comp_res++;
				}
		    		while (left <= right) {
		        		tmp_sum += join_size[left];
		        		if (tmp_sum > threshold){
		            			break;
					}
		        		map_idx = join_index[left];
		        		pos_con_map[map_idx] = free_map;
		        		left++;
		    		}
		    		free_map++;
			}
		}

    	}
	//if(flag){
	//	num_comp_res++;
	//}
	
	opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp_res*sizeof(opt_pk_t *));
	unsigned short int k1 = 0;
	array_comp_list_t * res = create_array_comp_list();
	cl = acl->head;
        size_t nbcomp = num_comp_res;
	
	unsigned short int ind = num_comp_res;
	for(k=0; k < num_comp; k++){
		unsigned short int k2 = num_comp_res - k1 - 1;
		
		unsigned short int comp_size = cl->size;
		if(pos_con_map[k]==2){
			poly[k2] = opt_poly_alloc(comp_size,0);
			poly[k2]->C = poly2[k]->C;
			poly[k2]->F = poly2[k]->F;
			opt_matrix_sort_rows(opk,poly[k2]->F);
			poly[k2]->nbeq = poly2[k]->nbeq;
			poly[k2]->nbline = poly2[k]->nbline;
			//poly[num_comp_res - k1 - 1]->satC = poly2[k]->satC;
			poly[k2]->satF = opt_satmat_alloc(poly[k2]->C->nbrows,opt_bitindex_size(poly[k2]->F->nbrows));
			combine_satmat(opk,poly[k2],comp_size,poly[k2]->F->nbrows,false);
			insert_comp_list(res,copy_comp_list(cl));
			
			k1++;
		}
		else if(pos_con_map[k]==4) {
			
            		poly[k2] = opt_poly_alloc(comp_size, 0);
            		poly[k2]->C = opt_matrix_copy(poly1[k]->C);
            		poly[k2]->nbeq = poly1[k]->nbeq;
            		size_t begin = poly1[k]->F->nbrows - nblinemapa[k];
            		opt_matrix_t *joinedF = opt_matrix_append(poly1[k]->F, poly2[k]->F);
            		poly[k2]->F = joinedF;
            		// remove_common_gen(opk,poly[k2]->F,begin);
            		poly[k2]->nbline = poly1[k]->nbline + poly2[k]->nbline;
            		opt_matrix_sort_rows(opk, poly[k2]->C);
            		size_t num = poly[k2]->F->nbrows - begin;
            		poly[k2]->satF = opt_satmat_alloc(poly[k2]->C->nbrows, opt_bitindex_size(poly[k2]->F->nbrows));
            		combine_satmat(opk, poly[k2], comp_size, begin, false);
            		opt_poly_dual(poly[k2]);
            		opt_cherni_add_and_minimize(opk, false, poly[k2], begin);
            		opt_poly_dual(poly[k2]);
			
            		/* max flow for deleting constraints */
			if (!opk->exn && (poly[k2]->C->nbrows > 1)) {
				array_comp_list_t * facl = compute_finest_partition(poly[k2]->C,cl,maxcols, NULL);
				
				unsigned short int *ca = to_sorted_array(cl,maxcols);
				comp_list_t *fcl = facl->head;
				
				nbcomp+= facl->size;
				
				poly = (opt_pk_t **)realloc(poly,(nbcomp-1)*sizeof(opt_pk_t*));
				
				opt_pk_t * tmp = poly[k2];
				opt_matrix_t * F = tmp->F;
				opt_matrix_t * C = tmp->C;
				
				bool first = false;	
				while(fcl!=NULL){
					unsigned short int * ca1 = to_sorted_array(fcl,maxcols);
					unsigned short int * ind_map_a = map_index(ca1,ca,fcl->size);
					comp_list_t * newCl = copy_comp_list(fcl);
					opt_matrix_t * C1 =  opt_matrix_alloc(C->nbrows+1,fcl->size+opk->dec,false);
					
					size_t nbeq = split_matrix(opk,C1,C,ind_map_a,fcl->size);
					
					unsigned short int ind1=0;					
					if(fcl->size<=threshold){
						
						if(!first){
							ind1 = k2;
						}
						else{
							ind1 = ind;
						}
						poly[ind1] = opt_poly_alloc(fcl->size,0);
						poly[ind1]->C = C1;
						poly[ind1]->nbeq =nbeq;
						poly[ind1]->F = opt_matrix_alloc(F->nbrows,fcl->size+opk->dec,false);
						poly[ind1]->nbline = split_matrix(opk,poly[ind1]->F,F,ind_map_a,fcl->size); 
						
						poly[ind1]->satC = opt_satmat_alloc(poly[ind1]->F->nbrows,opt_bitindex_size(poly[ind1]->C->nbrows));
						combine_satmat(opk,poly[ind1],fcl->size,poly[ind1]->C->nbrows,true);	
						free(ca1);
						
						if(!first){
							insert_comp_list(res, newCl);
							first = true;
						}
						else{
							ind++;
							insert_comp_list_tail(res,newCl);
						}
					}
					else{
						
						//poly[ind1] = opt_poly_alloc(fcl->size,0);
						//poly[ind1]->C = C1;
						//printf("Needed here\n");
						//print_comp_list(newCl,maxcols);	
						//fflush(stdout);
						if(newCl->size>threshold){
			     	   			//unsigned short int *sortedVars = to_sorted_array(newCl, maxcols);
			     	   			//int *deletedConstraints = find_deleted_constraints(poly[ind1]->C->nbcolumns, (int) poly[ind1]->C->nbrows,
			                                  //                                 poly[ind1]->C->p);
							//size_t * deletedConstraints = constraints_to_delete(poly[ind1]->C,newCl,maxcols);
			     	   			//newCl = eliminate_constraint_rows(poly[ind1], poly[ind1]->C, deletedConstraints, sortedVars,
			                                  //                             opk);
					
							//free(sortedVars);
							//printf("mincut\n");
							//fflush(stdout);
							opt_matrix_t * newC = opt_matrix_alloc(C1->nbrows,C1->nbcolumns,false);
							//printf("maxflow input\n");
							//opt_matrix_fprint(stdout,C1);
							//fflush(stdout);
							array_comp_list_t * tacl = NULL;
							if(ss==UNARY_WEIGHT_DELETE){
								tacl = constraints_to_delete(C1,newC,newCl,maxcols, threshold);
							}
							else if(ss==STOER_MIN_CUT){
								tacl = constraints_to_delete_stoer_wagner(C1,newC,newCl,maxcols,threshold);
							}
							else{
								tacl = constraints_to_delete_binary(C1,newC,newCl,maxcols,threshold);
							}
							//printf("maxflow output\n");
							//opt_matrix_fprint(stdout,newC);
							//fflush(stdout);
							//print_array_comp_list(tacl,maxcols);
							nbcomp+=tacl->size;
							poly = (opt_pk_t **)realloc(poly,(nbcomp-1)*sizeof(opt_pk_t*));
							comp_list_t * tcl = tacl->head;
							while(tcl!=NULL){
								if(!first){
									ind1 = k2;
								}
								else{
									ind1 = ind;
								}
								
								unsigned short int * ca_t = to_sorted_array(tcl,maxcols);
								unsigned short int * ind_map_t = map_index(ca_t,ca1,tcl->size);
								comp_list_t * new_tcl = copy_comp_list(tcl);
								poly[ind1] = opt_poly_alloc(tcl->size,0);
								poly[ind1]->C = opt_matrix_alloc(newC->nbrows+1,tcl->size+opk->dec,false);
								split_matrix(opk,poly[ind1]->C,newC,ind_map_t,tcl->size);
								size_t nbrows = poly[ind1]->C->nbrows;
								poly[ind1]->C->p[nbrows][0] = 1;
								poly[ind1]->C->p[nbrows][1] = 1;
								poly[ind1]->C->nbrows++;
								
								poly[ind1]->F = NULL;
			     	   				poly[ind1]->satC = poly[ind1]->satF = NULL;
			           				opt_poly_chernikova(man, poly[ind1], "join-new");
								if(!first){
									insert_comp_list(res,new_tcl);
									first = true;
								}
								else{
									ind++;
									insert_comp_list_tail(res,new_tcl);
								}
								tcl = tcl->next;
							}
							opt_matrix_free(C1);
							opt_matrix_free(newC);
			     			}
						
						
					}
					
					
					fcl = fcl->next;	
				}
				free(ca);
				opt_poly_clear(tmp);
				
			} 
			else{
                		insert_comp_list(res, copy_comp_list(cl));
			 }
            		k1++;
		}
		cl = cl->next;
	}

	
	unsigned short int k2 = num_comp_res - k1 - 1;
	if(counter>=1){
		unsigned short int left = 0, right = counter - 1;
		if(ms==NO_MERGE){
			while(left<=right){
				size_t num_vertex1a = 0, num_vertex2b = 0;
				size_t nbF = 0, nbC = 0, nbeq = 0, nbline = 0, nblinea = 0;
				comp_list_t *clp = create_comp_list();
				char *clp_map = (char *) calloc(maxcols, sizeof(char));	
				unsigned short int kp = join_index[left];
		        	union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
				num_vertex1a = num_vertex1[kp];
				num_vertex2b = num_vertex2[kp];
				nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
				nbC = nbC + poly1[kp]->C->nbrows;
				nbeq = nbeq + poly1[kp]->nbeq;
				nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
				nblinea += nblinemapa[kp];
		        	create_joined_partition(man, res, poly, poly1, poly2, acl, clp, nbF, nbC, nbeq, nbline, num_vertex1a,
		                                num_vertex2b, nblinea, num_vertex1, num_vertex2, maxcols, nblinemapa, ca_arr, pos_con_map,
		                                opk, k2, pos_con_map[kp], nbcomp);
				if(opk->exn==ELINA_EXC_OVERFLOW){
					opk->exn = ELINA_EXC_NONE;
					nbcomp--;
				}
				left++;
				k2--;
				free(clp_map);
			}						
		}
		else{
			while (left <= right) {	
				size_t num_vertex1a = 0, num_vertex2b = 0;
				size_t nbF = 0, nbC = 0, nbeq = 0, nbline = 0, nblinea = 0;
				
				comp_list_t *clp = create_comp_list();
				char *clp_map = (char *) calloc(maxcols, sizeof(char));	
		    		if (left == right) {
		        		unsigned short int kp = join_index[left];
		        		union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
					
					if (!num_vertex1a) {
					    num_vertex1a = num_vertex1[kp];
					} else {
					    num_vertex1a *= num_vertex1[kp];
					}
					if (!num_vertex2b) {
					    num_vertex2b = num_vertex2[kp];
					} else {
					    num_vertex2b *= num_vertex2[kp];
					}
					nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
					nbC = nbC + poly1[kp]->C->nbrows;
					nbeq = nbeq + poly1[kp]->nbeq;
					nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
					nblinea += nblinemapa[kp];
		        		create_joined_partition(man, res, poly, poly1, poly2, acl, clp, nbF, nbC, nbeq, nbline, num_vertex1a,
		                                		num_vertex2b, nblinea, num_vertex1, num_vertex2, maxcols, nblinemapa, ca_arr, 
								pos_con_map, opk, k2, pos_con_map[kp], nbcomp);
					if(opk->exn==ELINA_EXC_OVERFLOW){
						opk->exn = ELINA_EXC_NONE;
						nbcomp--;
					}
			
		        		break;
		    		}
				unsigned short int sum_p = 0;
				unsigned short int tmp_sum = 0;
				unsigned short int kp =0;
				if(ms==MERGE_SMALL_LARGE){
		    			sum_p = join_size[right];
		    			tmp_sum = join_size[right];
		    			kp = join_index[right];
					union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
					if (!num_vertex1a) {
					num_vertex1a = num_vertex1[kp];
			    		} else {
						num_vertex1a *= num_vertex1[kp];
			    		}
			    		if (!num_vertex2b) {
						num_vertex2b = num_vertex2[kp];
			    		} else {
						num_vertex2b *= num_vertex2[kp];
			    		}
					nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
			
		    			nbC = nbC + poly1[kp]->C->nbrows;
		    			nbeq = nbeq + poly1[kp]->nbeq;
		    			nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
		    			nblinea += nblinemapa[kp];
					right--;
				}
		    		while (left <= right) {
		        		tmp_sum += join_size[left];
		        		if (tmp_sum > threshold) break;

		        		sum_p += join_size[left];
		        		kp = join_index[left];
		        		union_comp_list(clp, find_comp_by_index(acl->head, kp), clp_map);
					if (!num_vertex1a) {
					    num_vertex1a = num_vertex1[kp];
					} else {
					    num_vertex1a *= num_vertex1[kp];
					}
					if (!num_vertex2b) {
					    num_vertex2b = num_vertex2[kp];
					} else {
					    num_vertex2b *= num_vertex2[kp];
					}
					nbF = nbF + poly1[kp]->F->nbrows + poly2[kp]->F->nbrows;
				
					nbC = nbC + poly1[kp]->C->nbrows;
					nbeq = nbeq + poly1[kp]->nbeq;
					nbline = nbline + poly1[kp]->nbline + poly2[kp]->nbline;
					nblinea += nblinemapa[kp];
					left++;
		    		}
				
		    		create_joined_partition(man, res, poly, poly1, poly2, acl, clp, nbF, nbC, nbeq, nbline, num_vertex1a,
		                            num_vertex2b,
		                            nblinea, num_vertex1, num_vertex2, maxcols, nblinemapa, ca_arr, pos_con_map, opk,
		                            k2, pos_con_map[kp], nbcomp);
				if(opk->exn==ELINA_EXC_OVERFLOW){
					opk->exn = ELINA_EXC_NONE;
					nbcomp--;
				}
			
		    		k2--;
				free(clp_map);
			}
		}		
		/*if(opk->exn){
			opt_pk_t *tmp = poly[0];
			unsigned short int k1;
			for(k1=0; k1 < num_comp_res - 1; k1++){
				poly[k1] = poly[k1+1];
			}
			opt_poly_clear(tmp);
		}*/
	}
	
	
	if(destructive){
		for(k=0; k < num_compa; k++){
			opt_poly_clear(poly_a[k]);
		}
		array_comp_list_t * tmp = acla;
		free(poly_a);
		op->poly = poly;
		op->acl = res;
		free(tmp);
	}
	else{
		op->poly = poly;
		op->acl = res;
	}

	clock_t finished = clock();
	
	//cl = acl->head;
	//unsigned short int k2;
	//unsigned short int nc = num_comp;
	for(k=0; k< num_comp; k++){
		if(pos_con_map[k]==1){
			opt_poly_clear(poly2[k]);
		}
		else if(pos_con_map[k]==2){
			opt_poly_clear(poly1[k]);
		}
		else{
			opt_poly_clear(poly1[k]);
			opt_poly_clear(poly2[k]);
			
		}
		
		free(poly1[k]);
		free(poly2[k]);
		free(ca_arr[k]);
		free(var_map_a[k]);
		free(var_map_b[k]);
	}
	
	/*cl = op->acl->head;	
        for(k=0; k < op->acl->size; k++){
		if(cl->size>threshold && pos_con_map[k]!=2){
			printf("BREAK %d %u %d\n",cl->size,pos_con_map[k],k);
			print_comp_list(cl,maxcols);
			break;
		}
		cl = cl->next;
	}*/
	free(rmapa);
	free(rmapb);
	free(counterFa);
	free(counterFb);
	free(counterCa);
	free(counterCb);
	free(nbmapCa);
	free(nbmapCb);
	free(nbmapFa);
	free(nbmapFb);
	free(nbeqmapa);
	free(nbeqmapb);
	free(nblinemapa);
	free(nblinemapb);
	free(var_map_a);
	free(var_map_b);
	free(ca_arr);
	free(poly1);
	free(poly2);
	free(num_vertex1);
	free(num_vertex2);
	free(num_vertex_a);
	free(num_vertex_b);
	free(pos_con_map);
	free_array_comp_list(acl);
	free(join_index);
	free(join_size);
	//opt_matrix_fprint(stdout,poly[0]->C);
	//opt_matrix_fprint(stdout,poly[0]->F);
	//FILE *fp;
	//fp = fopen("/tmp/seahorn1.txt","a");
	
	
   //fprintf(fp,"%d %d\n",max,max_finest);
	//print_array_comp_list(op->acl,op->maxcols);
	//print_array_comp_list(acl_finest,op->maxcols);
   //fflush(stdout);
   //fclose(fp);
  // free_array_comp_list(acl_finest);
	
	
	//printf("After JOIN OUTPUT %d %d %d %d\n",opk->exn, destructive,nbcomp,num_comp_res);
	//print_array_comp_list(op->acl,op->maxcols);
	//elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,op);
        //elina_lincons0_array_fprint(stdout,&arr3,NULL);
	//elina_lincons0_array_clear(&arr3);
	
	//for(k=0;k < op->acl->size; k++ ){
	//	printf("pointer: %p %p\n",op->poly[k],op->poly[k]->C);
		//opt_matrix_fprint(stdout,op->poly[k]->C);
	//}
	//if(op->poly[k] && op->poly[k]->C){
	//	opt_matrix_fprint(stdout,op->poly[k]->C);
	//}
	//print_array_comp_list(op->acl,maxcols);
	//fflush(stdout);
	opk->exn = ELINA_EXC_NONE;
	elina_interval_t ** interval = opt_pk_to_box(man,op);
	//FILE *fp;
	//fp = fopen("/tmp/seahorn.txt","a");
	
	opt_numint_t reward = 0;
	num_semi_inf=0;
	num_exact = 0;
	num_bounded = 0;
	for(k=0; k < maxcols-2; k++){
		if(reward_map[k]){
		
			//elina_interval_print(interval[k]);
			if(!elina_interval_is_top(interval[k])){
				if(elina_scalar_infty(interval[k]->sup) || elina_scalar_infty(interval[k]->inf)){
						num_semi_inf++;
				}
				else{
					if(elina_scalar_equal(interval[k]->sup,interval[k]->inf)){
						num_exact++;
					}
					else{
						num_bounded++;
					}
				}
			}
		}
		elina_interval_free(interval[k]);
	}

	if(max_block_gt_0){
		int cpuCylces = (int) finished - started;

		reward = 3*num_exact + 2*num_bounded + num_semi_inf;


		PyObject* pFunc = PyObject_GetAttrString(pModule, "learn");


		PyObject* pArgs = PyTuple_New(2);

		PyObject* pValue = PyLong_FromLong(reward);
		PyTuple_SetItem(pArgs, 0, pValue);

		pValue = PyLong_FromLong(cpuCylces);
		PyTuple_SetItem(pArgs, 1, pValue);
		
		// PyObject_CallObject(pFunc, pArgs);

		// Py_END_ALLOW_THREADS

		PyObject_CallObject(pFunc, pArgs);

		PyGILState_Release(gstate);
	}

	
	

	

	

	// PyEval_SaveThread();

	// PyEval_ReleaseLock();



	if(size_a&&!destructive){
		opt_pk_free(man,tmp_a);
	}
	if(size_b){
		opt_pk_free(man,tmp_b);
	}
	
	
	return op;
}


typedef struct precision_t {
  unsigned short int top;
  unsigned short int half_top;
  unsigned short int points;
  double volume;
} precision_t;

int compare_precision(precision_t *p1, precision_t *p2){
	if(p1->points>p2->points){
		return 1;
	}
	else if(p1->points < p2->points){
		return 0;
	}
	else{
		if(p1->volume > p2->volume ){
			return 1;
		}
		else if(p1->volume < p2->volume){
			return 0;
		}
		else{
			if(p1->half_top >= p2->half_top){
				return 1;
			}
			else{
				return 0;
			}
		}
	}
}

void copy_precision(precision_t *p1, precision_t *p2){
	p1->top = p2->top;
	p1->half_top = p2->half_top;
	p1->points = p2->points;
	p1->volume = p2->volume;
}

void merge_precision(precision_t *arr, size_t *arr2, size_t l, size_t m, size_t r) {
    size_t i, j, k;
    size_t n1 = m - l + 1;
    size_t n2 = r - m;

    /* create temp arrays */
    size_t L2[n1], R2[n2];
    precision_t * L = (precision_t *)malloc(n1*sizeof(precision_t));
    precision_t * R = (precision_t *)malloc(n2*sizeof(precision_t));
    
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) {
        copy_precision(&L[i],&arr[l + i]);
        L2[i] = arr2[l + i];
    }
    for (j = 0; j < n2; j++) {
        copy_precision(&R[j],&arr[m + 1 + j]);
        R2[j] = arr2[m + 1 + j];
    }

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (compare_precision(L+i,R+j)) {
            copy_precision(&arr[k],&L[i]);
            arr2[k] = L2[i];
            i++;
        } else {
            copy_precision(&arr[k], &R[j]);
            arr2[k] = R2[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1) {
        copy_precision(&arr[k] , &L[i]);
        arr2[k] = L2[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2) {
        copy_precision(&arr[k], &R[j]);
        arr2[k] = R2[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort_precision(precision_t *arr, size_t *arr2, size_t l, size_t r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        size_t m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort_precision(arr, arr2, l, m);
        mergeSort_precision(arr, arr2, m + 1, r);

        merge_precision(arr, arr2, l, m, r);
    }
}

opt_pk_array_t* opt_pk_join(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob)
{
  opt_pk_internal_t *opk = opt_pk_init_from_manager(man, ELINA_FUNID_JOIN);
  /*if(oa->maxcols<=step_size){
	bool write_flag = false;
	//printf("JOIN INPUT1 %d\n",oa->maxcols);
	//elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa);
  	//elina_lincons0_array_fprint(stdout,&arr1,NULL);
  	//elina_lincons0_array_clear(&arr1);
  	//elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,ob);
  	//elina_lincons0_array_fprint(stdout,&arr2,NULL);
  	//elina_lincons0_array_clear(&arr2);
	//fflush(stdout);
	#if defined (TIMING)
		start_timing();
  	#endif
	opt_pk_array_t * best = opt_poly_join_gen(man,oa,ob,destructive,oa->maxcols,0,0,false,&write_flag);
	#if defined (TIMING)
		record_timing(join_time);
  	#endif
	//printf("JOIN OUTPUT\n");
  	//elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,best);
  	//elina_lincons0_array_fprint(stdout,&arr3,NULL);
  	//elina_lincons0_array_clear(&arr3);
  	//fflush(stdout);
	
	return best;
  }*/
  
  //opt_pk_array_t ** res = (opt_pk_array_t **)malloc(num_thresholds*9*sizeof(opt_pk_array_t*));
  //opt_pk_array_t *best;
  #if defined (TIMING)
	start_timing();
  #endif
  //printf("JOIN INPUT %d\n",oa->maxcols);
  //elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa);
  //elina_lincons0_array_fprint(stdout,&arr1,NULL);
  //elina_lincons0_array_clear(&arr1);
  //elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,ob);
  //elina_lincons0_array_fprint(stdout,&arr2,NULL);
  //elina_lincons0_array_clear(&arr2);
  //fflush(stdout);
  //precision_t * precision = (precision_t*)calloc(num_thresholds*9,sizeof(precision_t));
  //unsigned short int l = 0;
  //unsigned short int best_threshold = 0;
  //int best_ss = 0, best_ms=0;
  //size_t *cost = (size_t *)calloc(num_thresholds*9,sizeof(size_t));
  //size_t *index_cost = (size_t *)malloc(num_thresholds*9*sizeof(size_t));
  //size_t *index_precision = (size_t *)malloc(num_thresholds*9*sizeof(size_t));
  
 
  /*for(unsigned short int threshold = step_size; threshold < oa->maxcols;threshold=threshold+step_size){
	for(int ss=STOER_MIN_CUT; ss <= BINARY_WEIGHT_DELETE; ss++){
		for(int ms = NO_MERGE; ms <= MERGE_SMALL_LARGE;ms++){
			
			res[l] = opt_poly_join_gen(man,oa,ob,destructive,threshold,ss,ms,first, &write_flag);
			if(first){
				first = false;
			}
			elina_interval_t **interval = opt_pk_to_box(man,res[l]);
			precision[l].top=0;
			precision[l].half_top = 0;
			precision[l].points = 0;	
			precision[l].volume = 0;
			for(unsigned short int i=0; i < oa->maxcols - 2; i++){
				if(elina_interval_is_top(interval[i])){
					precision[l].top++;
				}
				else if(elina_scalar_infty(interval[i]->sup) || elina_scalar_infty(interval[i]->inf)){
					precision[l].half_top++;
				}
				else if(elina_scalar_equal(interval[i]->sup,interval[i]->inf)){
					precision[l].points++;
				}
				else{
					double inf, sup;
					elina_double_set_scalar(&inf,interval[i]->inf,GMP_RNDU);
					elina_double_set_scalar(&sup,interval[i]->sup,GMP_RNDU);
					double width = sup-inf;
					if(!precision[l].volume){
						precision[l].volume = width;
					}
					else{
						precision[l].volume = precision[l].volume*width;
					}
				}
				elina_interval_free(interval[i]);
			}
			free(interval);
			opt_pk_t ** poly = res[l]->poly; 
			array_comp_list_t * acl = res[l]->acl;
			unsigned short int num_comp = acl->size;
			comp_list_t *cl = acl->head;
			for(unsigned short int k=0; k<num_comp; k++){
				opt_matrix_t * F = poly[k]->F;
				unsigned short int comp_size = cl->size;
				
				if(comp_size*F->nbrows > cost[l]){
					cost[l] = comp_size*F->nbrows;
				}
				
				cl = cl->next;	
			}
			
			index_cost[l] = l;
			index_precision[l] = l;
			l++;
		}
	}
  }
  
  mergeSort_precision(precision,index_precision,0,num_thresholds*9-1);
 
  mergeSort_size_t(cost,index_cost,0,num_thresholds*9-1);
  
  size_t min_index_sum = 18*num_thresholds;
  
  for(l=0; l < num_thresholds*9; l++){
	size_t index = index_precision[l];
	unsigned short int l1 = 0;
	
	for(; l1 < num_thresholds*9; l1++){
		if(index_cost[l1]==index){
			break;
		}
	}
	
	if(l+l1 < min_index_sum){
		min_index_sum = l+l1;
		
		best_threshold = index/9;
		int tmp = index%9;
		best_ss = tmp/3;
		best_ms = tmp%3;
		best = res[index];
	}
  }
  if(write_flag){
	  fprintf(fp,"%u %d %d ",(best_threshold+1)*step_size,best_ss,best_ms);
	  fprintf(fp,"\n");
	  fflush(fp);
  }
  for(unsigned short int i = 0; i < num_thresholds *9; i++){
		if(res[i]!=best){
			opt_pk_free(man,res[i]);
		}
  }*/
  opt_pk_array_t *res = opt_poly_join_gen(man,oa,ob,destructive);
  
  //free(precision);
  //free(cost);
  //free(index_cost);
  //free(index_precision);
  #if defined (TIMING)
	record_timing(join_time);
  #endif
  //printf("JOIN OUTPUT\n");
  //elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,res);
  //elina_lincons0_array_fprint(stdout,&arr3,NULL);
  //elina_lincons0_array_clear(&arr3);
  //fflush(stdout);
  return res;
}


void opt_poly_meet(bool meet,
	       bool lazy,
	       elina_manager_t* man,
	       opt_pk_array_t* op, opt_pk_array_t* oa, opt_pk_array_t* ob){
	
	opt_pk_internal_t *opk = opt_pk_init_from_manager(man, ELINA_FUNID_MEET);
	array_comp_list_t *acla = oa->acl;
	unsigned short int num_compa = acla->size;
	array_comp_list_t *aclb = ob->acl;
	unsigned short int num_compb = aclb->size;
	unsigned short int maxcols = oa->maxcols;
	/*************************
		Compute union of independent components
	*************************/
	array_comp_list_t *acl = union_array_comp_list(acla, aclb, maxcols);
	unsigned short int num_comp = acl->size;
	opt_pk_t **poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
	size_t * nbeqmap = (size_t *)calloc(num_comp,sizeof(size_t));
	/*****************************
		Factor A according to union
	*****************************/	
	opt_pk_t ** poly_a = oa->poly;	
	unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));
	size_t * nbmapa = (size_t *)calloc(num_comp,sizeof(size_t));
	comp_list_t * cla = acla->head;
	size_t i;
	unsigned short int k,j;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * oak = poly_a[k];	
		opt_matrix_t * ocak = oak->C;
		short int ind = is_comp_list_included(acl,cla,maxcols);
		rmapa[k] = ind;
		nbmapa[ind] = nbmapa[ind] + ocak->nbrows;
		nbeqmap[ind] = nbeqmap[ind] + oak->nbeq;
		cla = cla->next;
	}

	/*****************************
		Factor B according to union
	*****************************/		
	opt_pk_t ** poly_b = ob->poly;
	unsigned short int * rmapb = (unsigned short int *)calloc(num_compb, sizeof(unsigned short int));
	size_t * nbmapb = (size_t *)calloc(num_comp,sizeof(size_t));
	comp_list_t * clb = aclb->head;
	for(k = 0; k < num_compb; k++){
		opt_pk_t * obk = poly_b[k];
		opt_matrix_t * ocbk = obk->C;
		short int ind = is_comp_list_included(acl,clb,maxcols);
		rmapb[k] = ind;
		nbmapb[ind] = nbmapb[ind] + ocbk->nbrows;
		nbeqmap[ind] = nbeqmap[ind] + obk->nbeq;
		//opt_matrix_sort_rows(opk,ocbk);
		clb = clb->next;
	}

	comp_list_t * cl = acl->head;
	size_t * counterC = (size_t *)calloc(num_comp, sizeof(size_t));
	unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
	for(k = 0; k < num_comp; k++){
		unsigned short int comp_size = cl->size;
		poly[k] = opt_poly_alloc(comp_size,0);
		poly[k]->C = opt_matrix_alloc(nbmapa[k] + nbmapb[k]+1,comp_size+2,false);
		poly[k]->C->p[0][0] = 1;
		poly[k]->C->p[0][1] = 1;
		poly[k]->nbeq = nbeqmap[k];
		counterC[k] = 1;
		ca_arr[k] = to_sorted_array(cl,maxcols);
		cl = cl->next;
	}

	cla = acla->head;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		opt_matrix_t * src_mat;
		opt_matrix_t * dst_mat;
		size_t * counter;
		src_mat = src->C;
		dst_mat = dst->C;
		counter = counterC;
		unsigned short int comp_size = cla->size;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t i1 = counter[ind];
		counter[ind] = counter[ind] + nbconsa;
		for(i = 0; i < nbconsa; i++){
		    opt_numint_t * src_pi = src_p[i];
		    opt_numint_t * dst_pi = dst_p[i1];
		    dst_pi[0] = src_pi[0];
		    dst_pi[1] = src_pi[1];
		    unsigned short int l = 0;
		    for(j = 0; j < comp_size; j++){
			while(ca[l] != ca_a[j]){
			      l++;
			}
			dst_pi[l+2] = src_pi[j+2];
			l++;
		   }
		   i1++;
	        }
		free(ca_a);
		cla = cla->next;
	}
	
	clb = aclb->head;
	for(k = 0; k < num_compb; k++){
		opt_pk_t * src = poly_b[k];
		unsigned short int ind = rmapb[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_b = to_sorted_array(clb, maxcols);
		unsigned short int * ca = ca_arr[ind];
		opt_matrix_t * src_mat = src->C;
		opt_matrix_t * dst_mat = dst->C;
		unsigned short int comp_size = clb->size;
		size_t nbconsb = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t *counter = counterC;
		size_t i1 = counter[ind];
		counter[ind] = counter[ind] + nbconsb;
		for(i = 0; i < nbconsb; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[i1];
			dst_pi[0] = src_pi[0];
			dst_pi[1] = src_pi[1];
			unsigned short int l = 0;
			for(j = 0; j < comp_size; j++){
				while(ca[l] != ca_b[j]){
					l++;
				}
				dst_pi[l+2] = src_pi[j+2];
				l++;
			}
			i1++;
		}
		free(ca_b);
		clb = clb->next;
	}
	
	/*****************************
		Sort rows of each block of op
	*****************************/
	bool is_bottom = false;
	for(k=0; k < num_comp; k++){
		opt_pk_t * src = poly[k];
		opt_poly_chernikova(man,src,"meet abstract");
		if(src->F==NULL){
		   opt_poly_set_bottom(opk,op);
		   is_bottom = true;
		   break;
		}
	}
	if(!is_bottom){
		array_comp_list_t * tmp = oa->acl;
		if(op==oa){
			free_array_comp_list(tmp);
		}
		op->acl = acl;
		op->poly = poly;
	}
	else{
		for(k = 0; k < num_comp; k++){
			opt_pk_t * op_k = poly[k];
			if(op_k){
				opt_poly_clear(op_k);
			}
		}
		free(poly);
		free_array_comp_list(acl);
	}

	for(k=0; k < num_comp; k++){
		free(ca_arr[k]);
	}
	
	free(rmapa);
	free(rmapb);
	free(nbmapa);
	free(nbmapb);
	free(nbeqmap);
	free(ca_arr);
	free(counterC);
	//#if defined (CONVERT)
	//	free(nbgenmap);
	//	free(nblinemap);
	//	free(counterF);
	//	free(counterS);
	//	free(colmapS);
	//#endif
}


/* ********************************************************************** */
/* II. Meet */
/* ********************************************************************** */

/* ********************************************************************** */
/* II.1 Meet of two or more polyhedra */
/* ********************************************************************** */

opt_pk_array_t* opt_pk_meet_cons(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob){
	opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_MEET);
	opt_pk_array_t* op = destructive ? oa :  opt_pk_array_alloc(NULL,NULL,oa->maxcols);
	opt_pk_t ** poly_a = oa->poly;
	array_comp_list_t * acla = oa->acl;
	array_comp_list_t * aclb = ob->acl;

	if(oa->is_bottom || !acla){
		if(destructive){
			return oa;
		}
		else{
			free(op);
			return opt_pk_bottom(man,oa->maxcols - 2, 0);
		}
        }
	if(ob->is_bottom || !aclb){
		if(destructive){
			opt_poly_set_bottom(opk,oa);
			return oa;
		}
		else{
			free(op);
			return opt_pk_bottom(man,oa->maxcols - 2, 0);
		}
	}
	unsigned short int num_compa = acla->size;
	unsigned short int k;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * oak = poly_a[k];
		if(opk->funopt->algorithm>=0){
		   opt_poly_chernikova(man,oak,"meet abstract");
		}
		else{	
		    opt_poly_obtain_C(man,oak,"meet abstract");
		}
		if(opk->exn){
		   opk->exn = ELINA_EXC_NONE;
		   if(destructive){
			opt_poly_array_clear(opk,op);
			free(op);
		   }
		   else{
		   	free(op);
		   }
		   op = opt_pk_copy(man,ob);
		   return op;
		}
		if(!oak->C){
		   opt_poly_set_bottom(opk,op);
		   return op;
		}
	}
	
	opt_pk_t ** poly_b = ob->poly;
	
	unsigned short int num_compb = aclb->size;
	for(k = 0; k < num_compb; k++){
	    opt_pk_t * obk = poly_b[k];
	    if(opk->funopt->algorithm>=0){
	       opt_poly_chernikova(man,obk,"meet abstract");
	    }
	    else{	
		opt_poly_obtain_C(man,obk,"meet abstract");
	    }
	    if(opk->exn){
	       opk->exn = ELINA_EXC_NONE;
	       if(destructive){
		  return oa;
	       }
	       else{
		  free(op);
		  op = opt_pk_copy(man,oa);
		  return op;
	       }
	   }
	   if(!obk->C){
	      opt_poly_set_bottom(opk,op);
	      return op;
	   }
	}
	
	opt_poly_meet(true, opk->funopt->algorithm < 0,
		  man, op,oa,ob);
	
	return op;
}


opt_pk_array_t* opt_pk_meet(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob){
	#if defined (TIMING)
		start_timing();
	#endif
	opt_pk_array_t *op = opt_pk_meet_cons(man,destructive,oa,ob);
	#if defined (TIMING)
		record_timing(meet_time);
	#endif
	return op;
}




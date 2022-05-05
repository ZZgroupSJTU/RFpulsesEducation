// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eOpCube
//! @{



template<typename T1, typename eop_type>
class eOpCube : public BaseCube<typename T1::elem_type, eOpCube<T1, eop_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool prefer_at_accessor = ProxyCube<T1>::prefer_at_accessor;
  
  arma_aligned const ProxyCube<T1> P;
  
  arma_aligned const elem_type aux;        //!< storage of auxiliary data, user defined format
  arma_aligned const u32       aux_u32_a;  //!< storage of auxiliary data, u32 format
  arma_aligned const u32       aux_u32_b;  //!< storage of auxiliary data, u32 format
  arma_aligned const u32       aux_u32_c;  //!< storage of auxiliary data, u32 format
  
  inline         ~eOpCube();
  inline explicit eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          eOpCube(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices);
  
  arma_inline u32 get_n_rows()       const;
  arma_inline u32 get_n_cols()       const;
  arma_inline u32 get_n_elem_slice() const;
  arma_inline u32 get_n_slices()     const;
  arma_inline u32 get_n_elem()       const;
  
  arma_inline elem_type operator[] (const u32 i)                                   const;
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const;
  };



//! @}

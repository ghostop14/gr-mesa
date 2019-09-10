/* -*- c++ -*- */
/* 
 * Copyright 2019 ghostop14.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_MESA_IOSELECTOR_IMPL_H
#define INCLUDED_MESA_IOSELECTOR_IMPL_H

#include <mesa/ioselector.h>

namespace gr {
  namespace mesa {

    class ioselector_impl : public ioselector
    {
     private:
      // Nothing to declare in this block.
    	int d_numinputs;
    	int d_numoutputs;
    	int d_curInput;
    	int d_curOutput;
    	int d_itemsize;

     public:
      ioselector_impl(int numinputs,int numoutputs,int inputport,int outputport,int itemsize);
      ~ioselector_impl();

      virtual void set_input_index(int newValue);
      virtual void set_output_index(int newValue);

      void handleMsgInputIndex(pmt::pmt_t msg);
      void handleMsgOutputIndex(pmt::pmt_t msg);

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_IOSELECTOR_IMPL_H */


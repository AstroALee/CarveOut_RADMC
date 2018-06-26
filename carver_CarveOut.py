"""
   carver_CarveOut.py

   Purpose:
        Carving class that pulls the relevant layers needed to create a RadMC amr
        and number density files. Should not need to edit this file ever.

   Author:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.6.4, YT 3.4.1
"""

import numpy as np
import inputs_CarveOut as inputs
import yt
from yt.utilities.lib.write_array import write_3D_array, write_3D_vector_array


class RadMC3DLayer:
    '''

    This class represents an AMR "layer" of the style described in
    the radmc3d manual. Unlike yt grids, layers may not have more
    than one parent, so level L grids will need to be split up
    if they straddle two or more level L - 1 grids.

    '''
    def __init__(self, level, parent, unique_id, LE, RE, dim, is_periodic):
        self.level = level
        self.parent = parent
        self.LeftEdge = LE
        self.RightEdge = RE
        self.ActiveDimensions = dim
        self.id = unique_id
        self.is_periodic = is_periodic


    def get_overlap_with(self, grid):
        '''

        Returns the overlapping region between two Layers,
        or a layer and a grid. RE < LE means in any direction
        means no overlap.

        '''
        LE = np.maximum(self.LeftEdge,  grid.LeftEdge)
        RE = np.minimum(self.RightEdge, grid.RightEdge)
        return LE, RE

    def overlaps(self, grid):
        '''

        Returns whether or not this layer overlaps a given grid

        '''
        LE, RE = self.get_overlap_with(grid)
        if np.any(RE <= LE):
            return False
        else:
            return True



class CarvingWriter:
    '''
        Parameters:
        a_ds = pointer to Orion2 data read in from YT
        a_box[Left,Right] = triple arrays for bottom-left and upper-right corners of carve region
        a_Domain[Left,Right] = triple arrays for bottom-left and upper-right of entire domain (needed for periodic sides)
        a_boxDim = triple array giving number of cells in each dimension for layer
        a_max_level = max AMR level we will use
        a_is_periodic = 0,1 whether we are allowing for periodic domains

    '''
    def __init__(self, a_ds, a_boxLeft, a_boxRight, a_DomainLeft, a_DomainRight, a_boxDim, a_max_level, a_is_periodic):
        self.max_level = a_max_level
        self.cell_count = 0
        self.layers = []
        self.domain_dimensions = a_boxDim #ds.domain_dimensions
        self.domain_left_edge  = yt.YTArray(a_boxLeft,'cm') # assumes you passed in CGS #ds.domain_left_edge
        self.domain_right_edge = yt.YTArray(a_boxRight,'cm') # assumes you passed in CGS #ds.domain_right_edge
        self.fulldomain_left = a_DomainLeft # pulled from yt, already a YT array
        self.fulldomain_right = a_DomainRight # pulled from yt, already a YT array
        self.grid_filename = inputs.out_afname #'amr_grid.inp'
        self.ds = a_ds
        self.is_periodic = a_is_periodic

        base_layer = RadMC3DLayer(0, None, 0,
                                  self.domain_left_edge,
                                  self.domain_right_edge,
                                  self.domain_dimensions,
                                  self.is_periodic)

        self.layers.append(base_layer)
        self.cell_count += np.product(a_boxDim) #np.product(a_ds.domain_dimensions)

        # Sort the grids by level
        # This is a list of pointers. Changing values in sorted grid also changes
        # things in self.ds.index.grids (equivalently a_ds.index.grids)
        sorted_grids = sorted(a_ds.index.grids, key=lambda x: x.Level)

        # Add relevant grids to a master list (ignoring levels above our max level)
        for grid in sorted_grids:
            if grid.Level <= self.max_level:
                self._add_grid_to_layers(grid)

        # Done with constructor

    def _get_parents(self, grid):
        parents = []
        for potential_parent in self.layers: # you've previously sorted by level, so the parent has to be in this smaller list
            if potential_parent.level == grid.Level - 1:
                if potential_parent.overlaps(grid): # overlaps include periodic shifts
                    parents.append(potential_parent)
        return parents

    def _add_grid_to_layers(self, grid):
        parents = self._get_parents(grid)
        for parent in parents:
            LE, RE = parent.get_overlap_with(grid)
            N = (RE - LE) / grid.dds
            N = np.array([int(n + 0.5) for n in N])
            new_layer = RadMC3DLayer(grid.Level, parent.id,
                                     len(self.layers),
                                     LE, RE, N, self.is_periodic)
            self.layers.append(new_layer)
            self.cell_count += np.product(N)

    def write_amr_grid(self):
        '''
        This routine writes the "amr_grid.inp" file that describes the mesh
        radmc3d will use.

        '''
        dims = self.domain_dimensions
        LE = self.domain_left_edge
        RE = self.domain_right_edge

        # Taken from YT, fairly certain shouldn't be necessary, since O2 code_length is CGS
        # RadMC-3D wants the cell wall positions in cgs. Convert here:
        LE_cgs = LE.in_units('cm').d  # don't write the units, though
        RE_cgs = RE.in_units('cm').d

        # calculate cell wall positions
        xs = [str(x) for x in np.linspace(LE_cgs[0], RE_cgs[0], dims[0]+1)]
        ys = [str(y) for y in np.linspace(LE_cgs[1], RE_cgs[1], dims[1]+1)]
        zs = [str(z) for z in np.linspace(LE_cgs[2], RE_cgs[2], dims[2]+1)]

        # writer file header
        grid_file = open(self.grid_filename, 'w')
        grid_file.write('1 \n')  # iformat is always 1
        if self.max_level == 0:
            grid_file.write('0 \n')
        else:
            grid_file.write('10 \n')  # only layer-style files are supported
        grid_file.write('1 \n')  # only cartesian coordinates are supported
        grid_file.write('0 \n')
        grid_file.write('{}    {}    {} \n'.format(1, 1, 1))  # assume 3D
        grid_file.write('{}    {}    {} \n'.format(dims[0], dims[1], dims[2]))
        if self.max_level != 0:
            s = str(self.max_level) + '    ' + str(len(self.layers)-1) + '\n'
            grid_file.write(s)

        # write base grid cell wall positions
        for x in xs:
            grid_file.write(x + '    ')
        grid_file.write('\n')

        for y in ys:
            grid_file.write(y + '    ')
        grid_file.write('\n')

        for z in zs:
            grid_file.write(z + '    ')
        grid_file.write('\n')

        # write information about fine layers, skipping the base layer:
        for layer in self.layers[1:]:
            p = layer.parent
            dds = (layer.RightEdge - layer.LeftEdge) / (layer.ActiveDimensions)
            if p == 0:
                ind = (layer.LeftEdge - LE) / (2.0*dds) + 1
            else:
                parent_LE = np.zeros(3)
                for potential_parent in self.layers:
                    if potential_parent.id == p:
                        parent_LE = potential_parent.LeftEdge
                ind = (layer.LeftEdge - parent_LE) / (2.0*dds) + 1
            ix = int(ind[0]+0.5)
            iy = int(ind[1]+0.5)
            iz = int(ind[2]+0.5)
            nx, ny, nz = layer.ActiveDimensions / 2
            s = '{}    {}    {}    {}    {}    {}    {} \n'
            #s = s.format(p, ix, iy, iz, nx, ny, nz)
            s = s.format(p, ix, iy, iz, int(nx), int(ny), int(nz))
            grid_file.write(s)

        grid_file.close()

    def _write_layer_data_to_file(self, fhandle, field, level, LE, dim):
        cg = self.ds.covering_grid(level, LE, dim, num_ghost_zones=1)
        if isinstance(field, list):
            data_x = cg[field[0]]
            data_y = cg[field[1]]
            data_z = cg[field[2]]
            write_3D_vector_array(data_x, data_y, data_z, fhandle)
        else:
            data = cg[field]
            write_3D_array(data, fhandle)

    def write_line_file(self, field, filename):
        '''
        This method writes out fields in the format radmc3d needs to compute
        line emission.

        Parameters
        ----------

        field : string or list of 3 strings
            If a string, the name of the field to be written out. If a list,
            three fields that will be written to the file as a vector quantity.
        filename : string
            The name of the file to write the data to. The filenames radmc3d
            expects for its various modes of operation are described in the
            radmc3d manual.

        '''
        fhandle = open(filename, 'w')

        # write header
        fhandle.write('1 \n')
        fhandle.write(str(self.cell_count) + ' \n')

        # now write layers:
        # accesses actual state data in _write_layer_data_to_file
        for layer in self.layers:
            lev = layer.level
            if lev == 0:
                LE = self.domain_left_edge
                N = self.domain_dimensions
            else:
                LE = layer.LeftEdge
                N = layer.ActiveDimensions

            self._write_layer_data_to_file(fhandle, field, lev, LE, N)

        fhandle.close()

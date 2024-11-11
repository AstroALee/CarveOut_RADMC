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
#import math

class RadMC3DLayer:
    '''

    This class represents an AMR "layer" of the style described in
    the radmc3d manual. Unlike yt grids, layers may not have more
    than one parent, so level L grids will need to be split up
    if they straddle two or more level L - 1 grids.

    '''
    def __init__(self, level, parent, unique_id, LE, RE, dim, is_periodic, a_patches):
        self.level = level
        self.parent = parent
        self.LeftEdge = LE    # should not be adjusted if periodic
        self.RightEdge = RE   # should not be adjusted if periodic
        self.ActiveDimensions = dim
        self.id = unique_id
        self.is_periodic = is_periodic # 1 if periodic adjustments were made
        self.patches = a_patches # if periodic adjustments


    def get_overlap_with(self, grid):
        '''

        Returns the overlapping region between two Layers,
        or a layer and a grid. RE < LE means in any direction
        means no overlap.

        '''
        LE = np.maximum(self.LeftEdge,  grid.LeftEdge)
        RE = np.minimum(self.RightEdge, grid.RightEdge)

        if(self.is_periodic and self.level==0):
            # adjust for potential periodic shifts
            for shifted_grid in self.patches:
                ledge_shift, redge_shift = shifted_grid
                LE = np.maximum(yt.YTArray(ledge_shift,'cm'),  grid.LeftEdge)
                RE = np.minimum(yt.YTArray(redge_shift,'cm'), grid.RightEdge)
                if np.any(RE > LE):
                    break # ASSUMPTION: it can overlap with only one of the broken up boxes
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
    def __init__(self, a_ds, a_boxLeft, a_boxRight, a_domainPatches, a_DomainLeft, a_DomainRight, a_boxDim, a_max_level, a_allow_periodic, a_force_unigrid):
        self.max_level = a_max_level
        self.cell_count = 0
        self.layers = []
        self.domain_dimensions = a_boxDim #ds.domain_dimensions
        self.domain_left_edge  = yt.YTArray(a_boxLeft,'cm') # assumes you passed in CGS #ds.domain_left_edge
        self.domain_right_edge = yt.YTArray(a_boxRight,'cm') # assumes you passed in CGS #ds.domain_right_edge
        self.fulldomain_left = a_DomainLeft # pulled from yt, already a YT array
        self.fulldomain_right = a_DomainRight # pulled from yt, already a YT array
        self.grid_filename = inputs.out_afname #'amr_grid.inp'
        self.ds = a_ds # yt dataset
        self.allow_periodic = a_allow_periodic
        self.force_unigrid = a_force_unigrid
        self.domainPatches = a_domainPatches # domain to carve out broken into patches


        if(self.force_unigrid):
            #self.max_level = 0
            newDim      = [(2**a_max_level)*x for x in a_boxDim]
            self.domain_dimensions = newDim
            base_layer = RadMC3DLayer(0, None, 0,
                                      self.domain_left_edge,
                                      self.domain_right_edge,
                                      newDim,
                                      self.allow_periodic, self.domainPatches)
            self.cell_count += np.product(newDim) #np.product(a_ds.domain_dimensions)
        else:
            base_layer = RadMC3DLayer(0, None, 0,
                                  self.domain_left_edge,
                                  self.domain_right_edge,
                                  self.domain_dimensions,
                                  self.allow_periodic, self.domainPatches)
            self.cell_count += np.product(a_boxDim) #np.product(a_ds.domain_dimensions)



        self.layers.append(base_layer)


        # Sort the grids by level
        # This is a list of pointers. Changing values in sorted grid also changes
        # things in self.ds.index.grids (equivalently a_ds.index.grids)
        # Outputs require sorting, so we sort it here
        sorted_grids = sorted(a_ds.index.grids, key=lambda x: x.Level)

        # Add relevant grids to a master list (ignoring levels above our max level)
        if(self.force_unigrid):
            print("Forcing unigrid") # This will skip adding additiona layers
        elif(inputs.force_nested):
            print("Sorting Layers so they're nested")
            self._created_nested_layers(sorted_grids)
        else:
            for grid in sorted_grids:
                if grid.Level <= self.max_level:
                    self._add_grid_to_layers(grid)

        # Statistics on layers included
        layerCounts = (self.max_level+1)*[0]
        for layer in self.layers:
            layerCounts[layer.level] = layerCounts[layer.level] + 1
        print("Number of layers for each level: " +  str(layerCounts))

        # If we are carving out a region with AMR, it is possible that some AMR
        # layers extend beyond the cutout region. If this is the case, the layer
        # was truncated to fit in the desired box. However, this is gonna screw
        # up the calculation of the cell size (length of layer / dimensions),
        # since the LE and RE values are calculated by only the overlap regions.
        # So here we loop over the layers and make sure they adhere to the right
        # cell size
        print("Original number of cells: " + str(self.cell_count))
        for layer in self.layers[1:]: # 0th entry is base layer
            # base and current level cell size
            base_cell = (self.domain_right_edge - self.domain_left_edge)/np.array(self.domain_dimensions)
            layer_cell = base_cell/pow(2.0,layer.level) # by construction of this loop, layer.level >= 1
            oldCellCount = np.product(layer.ActiveDimensions) # if we change the box, we have to change the total cell count

            if(self.layers.index(layer)==9):
                print("Prefix = " + str(layer.LeftEdge) + " " + str(layer.RightEdge) + " " + str(layer.ActiveDimensions))
                print(layer_cell)
            # Find the parent of the current layer
            parent_layer = None
            if(layer.is_periodic and layer.level==1):
                # periodic shift
                for shifted_grid in self.domainPatches:
                    ledge_shift, redge_shift = shifted_grid
                    LE = np.maximum(yt.YTArray(ledge_shift,'cm'),  layer.LeftEdge)
                    RE = np.minimum(yt.YTArray(redge_shift,'cm'), layer.RightEdge)
                    if np.any(RE > LE):
                        LE = yt.YTArray(ledge_shift,'cm')
                        RE = yt.YTArray(redge_shift,'cm')
                        N = (RE-LE)/base_cell
                        parent_layer = RadMC3DLayer(0, None,  # the 'layer' here is one of the disjoint patches
                                        0, LE, RE, N, self.allow_periodic, [])
                        break # ASSUMPTION: it can overlap with only one of the broken up boxes
            else:
                for potential_parent in self.layers:
                    if potential_parent.id == layer.parent:
                        parent_layer = potential_parent
                        break # if this list is huge, this could save a few tick tocks
            # if the layer was truncated, the layer's LE and RE will be exactly the same as the parent layer
            if( np.all(layer.LeftEdge > parent_layer.LeftEdge) and np.all(layer.RightEdge < parent_layer.RightEdge) ):
                # layer is fully inside and cell sizes will be correct
                continue
            else:
                # something was violated in the above
                if( np.any(layer.LeftEdge < parent_layer.LeftEdge)):
                    print("LEFT: I don't think this should ever happen...")
                    print(layer.LeftEdge , parent_layer.LeftEdge)
                    for i in range(len(layer.LeftEdge)):
                        layer.LeftEdge[i] = max(layer.LeftEdge[i] , parent_layer.LeftEdge[i] )
                    # replaces only elements where layer.LE < parent_layer.LE
                if( np.any(layer.RightEdge > parent_layer.RightEdge)):
                    print("RIGHT: I don't think this should ever happen...") # this actually will happen
                    #print(layer.RightEdge , parent_layer.RightEdge)
                    # replaces only elements where layer.LE < parent_layer.LE
                    for i in range(len(layer.RightEdge)):
                        #print(layer.LeftEdge[i] , layer.RightEdge[i])
                        #print(parent_layer.LeftEdge[i] , parent_layer.RightEdge[i])
                        layer.RightEdge[i] = min(layer.RightEdge[i] , parent_layer.RightEdge[i] )
                        assert(layer.RightEdge[i] > layer.LeftEdge[i])
                # however, truncation of the layer was done... cell counts and cell sizes likely wrong
                minRE = np.array( len(layer.RightEdge)*[0] )
                for i in range(len(layer.RightEdge)):
                    minRE[i] = min(layer.RightEdge[i],parent_layer.RightEdge[i])
                newDim = np.floor( (minRE - layer.LeftEdge.d)/layer_cell.d )
                # if newDim is odd, we throw away another cell to make it even
                for i in range(len(newDim)):
                    if(newDim[i]%2): #   = 1 if odd
                        #print(self.layers.index(layer))
                        newDim[i] = newDim[i] - 1
                        if(newDim[i]<2):
                            print("New dimension in " + str(i) + " direction < 2 : " + str(newDim[i]))
                            newDim[i] = 2 # this might f things up...
                        assert(newDim[i]>1) # not sure if newDim = 0 would cause issues # also, this assert doesn't seem to work?
                layer.RightEdge = layer.LeftEdge + ( newDim * layer_cell ) # 'layer' here is a pointer, the ultimate list will update
                layer.ActiveDimensions = np.array([int(x) for x in newDim])
                if(self.layers.index(layer)==9):
                    print("Postfix = " + str(layer.LeftEdge) + " "  + str(layer.RightEdge) + " " + str(layer.ActiveDimensions))
                #print(layer.LeftEdge , parent_layer.LeftEdge, newDim)
                for i in range(len(layer.ActiveDimensions)):
                    assert(layer.ActiveDimensions[i]>0)
                if(np.product(layer.ActiveDimensions)<0):
                    print("ISSUE: " + str(layer.ActiveDimensions))
                self.cell_count = self.cell_count - oldCellCount + np.product(layer.ActiveDimensions)
                #print("Adjusted layer from " + str(oldCellCount) + " to " + str(np.product(layer.ActiveDimensions)))
        print("New number of cells: " + str(self.cell_count))


        #print("Out of loop = " + str(self.layers[9].LeftEdge) + " " + str(self.layers[9].RightEdge) + " " + str(self.layers[9].ActiveDimensions))

        # Done with constructor

    def _get_parents(self, grid):
        parents = []
        for potential_parent in self.layers: # you've previously sorted by level, so the parent has to be in this smaller list
            if potential_parent.level == grid.Level - 1:
                if potential_parent.overlaps(grid): # overlaps here MUST account for periodic shifts, but only for grid.Level = 1
                    parents.append(potential_parent)
        return parents

    def _created_nested_layers(self,sorted_grids):
        # at this point you have only added the base layer (level = 0)
        for cur_lev in range(1,self.max_level+1): # for each level, find the single rectangle that would cover all the cur_level data
                                                  # each level then has exactly one entry in the layer list
                                                  # there can only be one parent then for each higher level
            # find the only possible parent
            pIDX = -1
            for parent in self.layers:
                if(parent.level == cur_lev-1):
                    pIDX = self.layers.index(parent)
                    #LEparent = parent.LeftEdge  #yt.YTArray([1e100,1e100,1e100],'cm')
                    #REparent = parent.RightEdge #yt.YTArray([-1e100,-1e100,-1e100],'cm')
                    pID   = parent.id
            dx = 0
            LEmin = yt.YTArray([1e100,1e100,1e100],'cm')
            REmax = yt.YTArray([-1e100,-1e100,-1e100],'cm')
            for grid in sorted_grids:
                if(grid.Level==cur_lev):
                    LE, RE = self.layers[pIDX].get_overlap_with(grid)
                    dx = grid.dds
                    N = (RE-LE)/dx
                    if(np.all(N>1)): # N > 0 means there is overlap
                        LEmin = np.minimum(LEmin,LE) # numpy does it entry by entry
                        REmax = np.maximum(REmax,RE)
            # At this point LEmin and REmax is the largest volume for grids that
            # also overlap with the parent. This is our layer
            N = (REmax-LEmin)/dx
            N = np.array([int(n + 0.5) for n in N])
            print("Cells on level " + str(cur_lev) + ": N = " + str(N))
            new_layer = RadMC3DLayer(cur_lev, pID,
                                 len(self.layers),
                                 LEmin, REmax, N, self.allow_periodic, [])
            self.layers.append(new_layer)
            self.cell_count += np.product(N)









    def _add_grid_to_layers(self, grid):
        parents = self._get_parents(grid) # looks in the already added layers (which includes/may only be the base layer)
        for parent in parents: # if potential parents were found, this is done, else skipped
            LE, RE = parent.get_overlap_with(grid)
            N = (RE - LE) / grid.dds
            if(np.any(N<2)):
                print("Skipping thin layer!  grid level = " + str(grid.Level))
            if(np.all(N>1)): # truncation and roundoff sometimes results in razor-thin layers
                N = np.array([int(n + 0.5) for n in N])
                new_layer = RadMC3DLayer(grid.Level, parent.id,
                                     len(self.layers),
                                     LE, RE, N, self.allow_periodic, [])
                self.layers.append(new_layer)
                self.cell_count += np.product(N)

    def write_amr_grid(self):
        '''
        This routine writes the "amr_grid.inp" file that describes the mesh
        radmc3d will use.

        '''
        dims = self.domain_dimensions # if force_unigrid, this should have been overwritten to the right value

        LE = self.domain_left_edge # carved out region
        RE = self.domain_right_edge

        #CellCount = np.product(self.domain_dimensions) # check

        # Taken from YT, fairly certain shouldn't be necessary, since O2 code_length is CGS
        # RadMC-3D wants the cell wall positions in cgs. Convert here:
        LE_cgs = LE.in_units('cm').d  # don't write the units, though
        RE_cgs = RE.in_units('cm').d  # also prevents passing by pointer (self.domain... etc. needs to be used again)

        # Shift so centered at 0,0,0
        Center_cgs = 0.5*(LE_cgs + RE_cgs)
        LE_cgs = LE_cgs - Center_cgs
        RE_cgs = RE_cgs - Center_cgs

        # calculate cell wall positions (may potentially exceed full domain if periodic)
        xs = [str(x) for x in np.linspace(LE_cgs[0], RE_cgs[0], dims[0]+1)]
        ys = [str(y) for y in np.linspace(LE_cgs[1], RE_cgs[1], dims[1]+1)]
        zs = [str(z) for z in np.linspace(LE_cgs[2], RE_cgs[2], dims[2]+1)]

        # writer file header
        grid_file = open(self.grid_filename, 'w')
        grid_file.write('1 \n')  # iformat is always 1
        if(self.max_level == 0 or self.force_unigrid==1):
            grid_file.write('0 \n')
        else:
            grid_file.write('10 \n')  # only layer-style files are supported
        grid_file.write('1 \n')  # only cartesian coordinates are supported
        grid_file.write('0 \n')
        grid_file.write('{}    {}    {} \n'.format(1, 1, 1))  # assume 3D
        grid_file.write('{}    {}    {} \n'.format(dims[0], dims[1], dims[2]))
        if(self.max_level != 0 and self.force_unigrid==0):
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

        # write information about fine layers
        for layer in self.layers[1:]: # [1:] skips entry 0, the base layer
            p = layer.parent # parent id
            dds = (layer.RightEdge - layer.LeftEdge) / (layer.ActiveDimensions) # cell size
            if p == 0: # parent is the base layer,
                # have to incorporate periodicity, if used
                # LE = base layer, needs to be adjusted
                # layer.Left is the level=1 layer, no adjustment needed
                # Beginning index
                ind = (layer.LeftEdge - LE) / (2.0*dds) + 1 # LE wasn't shifted so 0,0,0 is center. OK!
                if(layer.is_periodic):
                    for shifted_grid in self.domainPatches:
                        ledge_shift, redge_shift = shifted_grid
                        LE = np.maximum(yt.YTArray(ledge_shift,'cm'),  layer.LeftEdge)
                        RE = np.minimum(yt.YTArray(redge_shift,'cm'), layer.RightEdge)
                        if np.any(RE > LE):
                            print("HOW COULD THIS HAPPEN?!?!")
                            LE = yt.YTArray(ledge_shift,'cm')
                            ind = (layer.LeftEdge - LE) / (2.0*dds) + 1
                            #print(ind)
                            break # ASSUMPTION: it can overlap with only one of the broken up boxes
            else: # parent is an AMR layer
                parent_LE = np.zeros(3)
                for potential_parent in self.layers:
                    if potential_parent.id == p:
                        parent_LE = potential_parent.LeftEdge
                        break # might save a few tick tocks
                ind = (layer.LeftEdge - parent_LE) / (2.0*dds) + 1 # #index in parent grid cells; periodic or not, they both should be correct; difference is good here
            #print(str(np.array(ind)+0.5) + " " + str([int(x) for x in np.array(ind)+0.5]))
            ix = int(ind[0]+0.5) # beginning index in terms of parent grid cells (b/c of the 2*dds)
            iy = int(ind[1]+0.5)
            iz = int(ind[2]+0.5)
            #ix = math.ceil(ind[0]) # beginning index in terms of parent grid cells (b/c of the 2*dds)
            #iy = math.ceil(ind[1])
            #iz = math.ceil(ind[2])
            nx, ny, nz = layer.ActiveDimensions / 2 # number of cells (/2 to measure in dimensions of parent cell size)
            #print(nx,ny,nz)
            s = '{}    {}    {}    {}    {}    {}    {} \n'
            #s = s.format(p, ix, iy, iz, nx, ny, nz)
            #s = s.format(p, ix, iy, iz, int(np.rint(nx.d)), int(np.rint(ny.d)), int(np.rint(nz.d))) # RAD complains if anything is not int()
            s = s.format(p, ix, iy, iz, int(nx), int(ny), int(nz) ) # RAD complains if anything is not int()
            #s = s.format(p, ix, iy, iz, int(round(nx)), int(round(ny)), int(round(nz)) ) # RAD complains if anything is not int()
            #CellCount = CellCount + 2*(round(nx)*round(ny)*round(nz))
            grid_file.write(s)

        grid_file.close()
        #print("Cell count post AMR : " + str(CellCount))

    def _write_layer_data_to_file(self, fhandle, field, level, LE, dim):
        # coverting grids take into account periodicity, input original values for LE/RE

        lev = level
        if(self.force_unigrid):
            lev = self.max_level

        cg = self.ds.covering_grid(lev, LE, dim, num_ghost_zones=1)
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
                LE = self.domain_left_edge #passing to covering grid, periodicity incorporated on its own
                N = self.domain_dimensions
            else:
                LE = layer.LeftEdge        # passing to covering grid, pass original non-adjusted values
                N = layer.ActiveDimensions

            self._write_layer_data_to_file(fhandle, field, lev, LE, N)

        fhandle.close()

    def write_dust_file(self, field, filename):
        '''
        This method writes out fields in the format radmc3d needs to compute
        dust emission.

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
        fhandle.write('1 \n')

        # now write layers:
        # accesses actual state data in _write_layer_data_to_file
        print("Writing dust with self.cell_count = " + str(self.cell_count))
        CountUp = 0
        for layer in self.layers:
            lev = layer.level
            if lev == 0:
                LE = self.domain_left_edge #passing to covering grid, periodicity incorporated on its own
                N = self.domain_dimensions
            else:
                LE = layer.LeftEdge        # passing to covering grid, pass original non-adjusted values
                N = layer.ActiveDimensions

            self._write_layer_data_to_file(fhandle, field, lev, LE, N)
            CountUp = CountUp + np.product(N)
        print("Wrote dust with " + str(CountUp) + " (should match self.cell_count)")

        fhandle.close()

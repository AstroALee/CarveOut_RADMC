


import numpy as np
from inputs_CV import inputs
import yt
from yt.utilities.lib.write_array import write_3D_array, write_3D_vector_array


class RadMC3DLayer:
    '''

    This class represents an AMR "layer" of the style described in
    the radmc3d manual. Unlike yt grids, layers may not have more
    than one parent, so level L grids will need to be split up
    if they straddle two or more level L - 1 grids.

    '''
    def __init__(self, level, parent, unique_id, LE, RE, dim):
        self.level = level
        self.parent = parent
        self.LeftEdge = LE
        self.RightEdge = RE
        self.ActiveDimensions = dim
        self.id = unique_id


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
    def __init__(self, a_ds, a_boxLeftCell, a_DomainLeft, a_boxDim, a_boxDX):

        self.ds = a_ds # yt dataset
        # Entire domain
        self.domain_left = a_DomainLeft # left edge of domain, cgs
        self.domain_dim = a_boxDim # box dimensions, base grid
        self.domain_dx = yt.YTArray(a_boxDX,'cm') # box cell size, base grid, cgs

        # Carving region
        self.box_leftCell = a_boxLeftCell # carving region, lower-left corner grid cell dimension
        self.box_dim = inputs.Ncells # dimensions of carved out region in cells either on base grid or in forced unigrid
        self.max_level = inputs.max_level # AMR level to carve
        self.force_unigrid = inputs.is_unigrid

        # Initialized variables
        self.cell_count = 0
        self.layers = []
        self.grid_filename = inputs.out_afname #'amr_grid.inp'

        # If we are forcing unigrid, then the dimensions of the output is not the same as the base grid, potentially
        if(self.force_unigrid):
            newDim  = [(2**self.max_level)*x for x in self.box_dim]
            self.box_dim = newDim
            cellSize = self.domain_dx/(2**self.max_level)

            LE = self.domain_left + self.box_leftCell*self.domain_dx
            RE = LE + cellSize*self.box_dim

            # this is the entire layer for unigrid
            base_layer = RadMC3DLayer(0, None, 0,LE,RE,newDim)
            self.cell_count += np.product(newDim)

        # else if we not forcing unigrid, then we just use the base level
        else:
            LE = self.domain_left + self.box_leftCell*self.domain_dx
            RE = LE + self.domain_dx*self.box_dim
            # this is the base layer
            base_layer = RadMC3DLayer(0, None, 0,LE, RE, self.box_dim)
            self.cell_count += np.product(self.box_dim)

        self.layers.append(base_layer)

        # Sort the grids by level ( won't use this after layers are added to self)
        # This is a list of pointers. Changing values in sorted grid also changes
        # things in self.ds.index.grids (equivalently a_ds.index.grids)
        # Outputs require sorting, so we sort it here
        sorted_grids = sorted(self.ds.index.grids, key=lambda x: x.Level)
        print("length",len(sorted_grids))
        print("max level",self.max_level)

        # now looks through sorted grids and adds layers that are relevant
        # for what we want to carve out
        # Add relevant grids to a master list (ignoring levels above our max level)
        if(self.force_unigrid):
            print("Forcing unigrid") # This will skip adding additiona layers
        else:
            for grid in sorted_grids:
                if grid.Level <= self.max_level:
                    self._add_grid_to_layers(grid)

        # Statistics on layers included
        layerCounts = (self.max_level+1)*[0]
        for layer in self.layers:
            layerCounts[layer.level] = layerCounts[layer.level] + 1
        print("Number of layers for each level: " +  str(layerCounts))


    # end of constructor


    # children functions

    def _get_parents(self, grid):
        parents = []
        for potential_parent in self.layers: # you've previously sorted by level, so the parent has to be in this smaller list
            if potential_parent.level == grid.Level - 1:
                if potential_parent.overlaps(grid): # overlaps here MUST account for periodic shifts, but only for grid.Level = 1
                    parents.append(potential_parent)
        return parents

    def _add_grid_to_layers(self, grid):
        parents = self._get_parents(grid) # looks in the already added layers (which includes/may only be the base layer)
        for parent in parents: # if potential parents were found, this is done, else skipped
            LE, RE = parent.get_overlap_with(grid)
            N = (RE - LE) / grid.dds
            print("Add Layer LE/RE: ",LE.d,RE.d,N)
            if(np.any(N<2)):
                print("Skipping thin layer!  grid level = " + str(grid.Level))
            if(np.all(N>4)): # truncation sometimes results in razor-thin layers
                N = np.array([int(n + 0.5) for n in N])
                new_layer = RadMC3DLayer(grid.Level, parent.id,
                                     len(self.layers), LE, RE, N)
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

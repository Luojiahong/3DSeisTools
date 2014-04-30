"""
A submodule to provide core functionality and classes needed by various
software components in the 3DSeisTools package.
"""
if 'os' not in locals(): import os
if 'time' not in locals(): import time
from numpy import * #Should probably be more specific
#from matplotlib import pyplot as plt
#from array import array #This is for fast io when writing binary #AAA MAY NEED TO DELETE

#Read ASCII arrival time file output by FMM code

def num(s):
    """
    NEEDS TO BE UPDATED
    Convert a string to a number, choosing int or float automatically
    SLOW, don't use for large lists
     """
    try:
        return int(s)
    except ValueError:
        return float(s)

def parse_cfg(config_file):
    """
    Parse .cfg configuration file and return dictionary of contents.

    Arguments:
    config_file - Path to configuration file.

    Returns:
    mydict - Dictionary of parameters parsed from config_file.

    Example:
    In [1]: from anfseistools.core import parse_cfg

    In [2]: cfg_dict = parse_cfg('test_pf_2_cfg.cfg')

    In [3]: print cfg_dict
    {'misc': {'earth_radius': 6371.0,
              'tt_map_dir': '/Users/mcwhite/staging/tt_maps/June2010/'
             },
             'propagation_grid': {'minlon': -117.80,
                                  'dlon': 0.0327,
                                  'nlat': 76,
                                  'minlat': 32.5,
                                  'minz': 3.0,
                                  'dlat': 0.0273,
                                  'nr': 25,
                                  'dr': 2.0,
                                  'nlon': 73,
                                  'refinement_factor': 5,
                                  'ncells': 10
                                  },
             'location_parameters': {'buff1': 7,
                                     'buff2': 7,
                                     'dstep1': 5,
                                     'dstep2': 5,
                                     'nlat': 73,
                                     'nr': 25,
                                     'nlon': 76
                                     }
    }
    """
    import ConfigParser
    from antpy import eval_dict
    config = ConfigParser.RawConfigParser()
    config.read(config_file)
    mydict = {}
    for section in config.sections():
        section_dict = {}
        for option in config.options(section):
            section_dict[option] = config.get(section, option)
        mydict[section] = section_dict
    return eval_dict(mydict)

def load_faults():
    """
    NEEDS TO BE UPDATED
    Load the california fault map.
    """
    fnam = 'cal_faults.dat'
    fid = open(fnam,'r')
    a = fid.readlines()
    faultlon = [];faultlat = []
    for line in a:
        tmp = line.strip().split()
        if isnan(num(tmp[0])):
            faultlon.append(nan)
            faultlat.append(nan)
        else:
            faultlon.append(num(tmp[0]))
            faultlat.append(num(tmp[1]))
    print 'Faults loaded...'
    return faultlon, faultlat

def read_tt_vector(stanames, ind, ttdir):
    """
    NEEDS TO BE UPDATED
    Read from binary files to create a vector of traveltimes for the
    stations in stanames at the index location ind, which is the 1D
    index stanames is a list of station names only.
    """
    from numpy import array
    ttvec = array([])
    for sta in stanames:
        #fid = open(ttdir+'bin.'+sta+'.traveltime')
        file_path = '%sbin.%s.traveltime' % (ttdir, sta)
        if not os.path.isfile(file_path):
            return None
        fid = open('%sbin.%s.traveltime' % (ttdir, sta))
        ttvec = append(ttvec, read_binary_float(fid,ind) )
        fid.close()
    return ttvec

class Fmm_vgrids():
    """
    NEEDS TO BE UPDATED
    vgrids.in is the FMM velocity file format. Looks like:
    ngrids ntypes
    nradius nlat nlon     !number of points
    dradius dlat dlon     !grid spacings; uniform in each direction
    oradius olat olon     !origin of the grid
    V(1,1,1)
    V(2,1,1)
    V(1,2,1)
    V(1,1,2)
    etc...

    NOTE: Right-handed; oradius is somewhere deep in the earth

    ngrids is the number of volumes. We generally use 1
    ntypes will be 1 for just P, 2 for P and S
    """
    def __init__(self):
        """
        NEEDS TO BE UPDATED
        """
        pass

    def read(self,fnam='vgrids.in'):
        """
        NEEDS TO BE UPDATED
        Create a matrix of velocities.
        """
        self.fnam = fnam
        fid = open(fnam,'r')
        tmp = fid.readline().strip().split()
        self.ngrids,self.ntypes = num(tmp[0]),num(tmp[1])
        tmp = fid.readline().strip().split()
        self.nrad,self.nlat,self.nlon = num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp = fid.readline().strip().split()
        self.drad,self.dlat,self.dlon = num(tmp[0]),num(tmp[1]),num(tmp[2])
        tmp = fid.readline().strip().split()
        self.orad,self.olat,self.olon = num(tmp[0]),num(tmp[1]),num(tmp[2])
        #Loop to create a numpy matrix
        self.vel = empty((self.nlon,self.nlat,self.nrad))
        for irad in range(self.nrad):
            for ilat in range(self.nlat):
                for ilon in range(self.nlon):
                    tmp = fid.readline().strip().split()
                    self.vel[ilon,ilat,irad] = float(tmp[0])
        #THERE SHOULD BE MORE STATEMENTS HERE IN CASE NTYPES,NGRIDS!= 1
        fid.close()
    def write(self,outfnam='out.vgrids'):
        """
        NEEDS TO BE UPDATED
        Write to vgrids format.
        """
        self.outfnam = outfnam
        fid = open(outfnam,'w')
        outs = str(self.ngrids)+' '+str(self.ntypes)+'\n'
        fid.write(outs)
        outs = str(self.nrad)+' '+str(self.nlat)+' '+str(self.nlon)+'\n'
        fid.write(outs)
        outs = str(self.drad)+' '+str(self.dlat)+' '+str(self.dlon)+'\n'
        fid.write(outs)
        outs = str(self.orad)+' '+str(self.olat)+' '+str(self.olon)+'\n'
        fid.write(outs)
        for ilon in range(self.nlon):
            for ilat in range(self.nlat):
                for irad in range(self.nrad):
                    fid.write('{1:.{0}f}'.format(3,self.vel[ilon,ilat,irad])+'\n')
        fid.close()

class Interface():
    """
    NEEDS TO BE UPDATED
    A 2D matrix containing interface information, e.g., topography or
    moho depth.
    """
    def __init__(self):
        pass

    def set_data(self, x, y, z):
        """
        NEEDS TO BE UPDATED
        This builds the class ad-hoc. x,y,z are all numpy arrays. x,y
        are vectors, z is 2D.
        """
        self.x, self.y, self.z = x, y, z
        deg_to_rad = math.acos(-1) / 180.0  #math.acos(-1) is pi
        self.nx, self.ny = len(x), len(y) #total number of x and y coordinates
        self.ox, self.oy = min(x), min(y)   #origins in lat,lon
        self.dx = self.x[1] - self.x[0]     #Spacings in degrees
        self.dy = self.y[1] - self.y[0]
        self.ox_rad = self.ox * deg_to_rad  #Origins in radians
        self.oy_rad = self.oy * deg_to_rad
        self.dx_rad = self.dx * deg_to_rad  #Spacings in radians
        self.dy_rad = self.dy * deg_to_rad

    def read_topo_netcdf(self, fnam, subsamp=1):
        """
        NEEDS TO BE UPDATED
        Read a netcdf formatted topography file
        """
        from scipy.io import netcdf_file as netcdf
        topo = netcdf(fnam,'r')
        x = topo.variables['x'][::subsamp] #vector
        y = topo.variables['y'][::subsamp] #vector
        z = topo.variables['z'][::subsamp, ::subsamp] #matrix
        self.subsamp = subsamp
        self.set_data(x, y, z)

    def interp(self, xi, yi):
        """
        NEEDS TO BE UPDATED
        Find the values at points xi,yi by cubic spline
        """
        import scipy.interpolate
        gnn = scipy.interpolate.RectBivariateSpline(self.x,
                                                    self.y,
                                                    self.z.transpose())
        self.x_unint = self.x #Remember non-interpolated values
        self.y_unint = self.y
        self.z_unint = self.z
        zi = gnn.__call__(xi, yi) #the actual interpolation
        self.set_data(xi, yi, zi)

    def write_fmm(self, fnam='out_interfaces', ifappend=False):
        """
        NEEDS TO BE UPDATED
        Write to the fmm ascii format. ifappend should be set to True
        for any interface after the first; this will only write the z
        values and not the header
        """
        if not ifappend:
            fid = open(fnam, 'w')
            fid.write('2\n') # THIS ASSUMES THAT THERE ARE ONLY 2 INTERFACES
            fid.write('%u %u\n'%(self.ny,self.nx) )
            fid.write('%f %f\n'%(self.dy_rad,self.dx_rad) )
            fid.write('%f %f\n'%(self.oy_rad,self.ox_rad) )
        else:
            fid = open(fnam, 'a')
        for iy in range(self.ny):#Write the interface as a vector
            outs = ''.join('{1:.{0}f}\n'.format(3,ii)  for ii in self.z[:, iy])
            fid.write(outs)
        fid.close()

def find_containing_cube(px, py, pz, xvec, yvec, zvec):
    """
    NEEDS TO BE UPDATED
    Find the 8 endpoints for the cell which contains point px,py
    We take advantage of the regular grid
    Assumes the point is inside the volume defined by xvec,yvec,zvec
    Returns an array of size 8,3 where the rows contain x,y,z
    coordinates of the cubes endpoints
    Also returns indexes of endpoints
    """
    #Find the nearest node point and indexes <--"indices" sounds stupid to me
    xind, xnode = _find_nearest(px, xvec)
    yind, ynode = _find_nearest(py, yvec)
    zind, znode = _find_nearest(pz, zvec)
    #Now check if the 3 coordinates of p are greater or less than the
    #node it is nearest.
    if px >= xnode:
    #px is east of the nearest node
        #if px is on the x boundary, return a duplicate point
        xi, xn = (xind, xvec[xind]) if px == max(xvec)\
                else (xind + 1, xvec[xind + 1])
    else:
    #px is west of the nearest node
        #if px is on the x boundary, return a duplicate point
        xi, xn = (xind, xvec[xind]) if px == min(xvec)\
                else (xind - 1, xvec[xind - 1])
    if py >= ynode:
    #py is north of the nearest node
        #if py is on the y boundary, return a duplicate point
        yi, yn = (yind, yvec[yind]) if py == max(yvec)\
                else (yind + 1, yvec[yind + 1])
    else:
    #px is south of the nearest node
        yi, yn = (yind, yvec[yind]) if py == min(yvec)\
                else (yind - 1, yvec[yind - 1])
    if pz <= znode:
    #pz is above the nearest node
        zi, zn = (zind, zvec[zind]) if pz == max(zvec)\
                else (zind + 1, zvec[zind + 1])
    else:
    #pz is below the nearest node
        zi, zn = (zind, zvec[zind]) if pz == min(zvec)\
                else (zind - 1, zvec[zind - 1])
    #Add new endpoints to define the cube
    endpoints = []
    endpoints.append([xnode, ynode, znode])
    endpoints.append([xn, ynode, znode])
    endpoints.append([xn, yn, znode])
    endpoints.append([xnode, yn, znode])
    endpoints.append([xnode, ynode, zn])
    endpoints.append([xn, ynode, zn])
    endpoints.append([xn, yn, zn])
    endpoints.append([xnode, yn, zn])
    #Add indices
    indexes = []
    indexes.append([xind, yind, zind])
    indexes.append([xi, yind, zind])
    indexes.append([xi, yi, zind])
    indexes.append([xind, yi, zind])
    indexes.append([xind, yind, zi])
    indexes.append([xi, yind, zi])
    indexes.append([xi, yi, zi])
    indexes.append([xind, yi, zi])
    return endpoints, indexes

def find_nearest(nparray,value):
    """
    NEEDS TO BE UPDATED
    Returns the nearest item in nparray to value
    """
    idx = (abs(nparray-value)).argmin()
    return nparray.flat[idx]

def find_nearest_index(nparray,value):
    """
    NEEDS TO BE UPDATED
    """
    idx = (abs(nparray-value)).argmin()
    return idx

def _find_nearest(px, xvec):
    """
    NEEDS TO BE UPDATED
    Find the nearest x in xvec
    returns index
    """
    best_ind = 0
    shortest = 100000000.0;
    for ii in range(len(xvec)):
        if abs(xvec[ii]-px)<shortest:
            shortest = abs(xvec[ii]-px)
            best_ind = ii
    return best_ind, xvec[best_ind]

def read_binary_float(fid, n=0, precision='double'):
    """
    NEEDS TO BE UPDATED
    Read the nth float value from a binary file at the with the given
    precision following python conventions, 0 is the index of the first
    value.
    """
    import struct
    if precision is 'single':
        numbytes = 4;packstr = 'f'
    else:
        numbytes = 8;packstr = 'd'
    offset = n*numbytes #Go to the right spot in the file
    fid.seek(offset)
    raw = fid.read(numbytes)
    tmp = struct.unpack(packstr,raw)
    val = tmp[0]
    return val

def _grid_search_traveltimes_origin(arrsta, qx, qy, qz, arrvec, li):
    """
    NEEDS TO BE UPDATED
    Find the minimum value of the origin time standard deviation
    following Ben-Zion et al., 1992 (JGR)

    Arguments:
    arrsta - list of station names; strings
    qx, qy, qz - vectors of indices to search through
    arrvec - vector of absolute arrivals in the same order as sta
    li - LinearIndex class for the entire traveltime grid
    """
    #There is no reason this should be here, but the next line
    #generated error messages if it wasn't. I'm confused.
    from numpy import array
    rms = array([])
    #origin_std = array([])
    #Give large starting values
    origin_std = empty([ len(qy),len(qx),len(qz)] )+1000
    origin_mean = empty([ len(qy),len(qx),len(qz)] )
    search_inds = LinearIndex(len(qx),len(qy),len(qz))
    #Loop over the three vectors, searching every point
    for ix in range(len(qx)):
        for iy in range(len(qy)):
            for iz in range(len(qz)):
                #initialize the calculated tt vector
                calctt = array([])
                #Find the vector index
                ind = li.get_1D(qx[ix],qy[iy],qz[iz])
                for sta in arrsta: #Build vector of calculated ttimes
                    fid = open('bin.'+sta+'.traveltime')
                    tmp = read_binary_float(fid,ind)
                    #if tmp<0: tmp = NaN #Kill negative travel times
                    calctt = append(calctt, tmp )
                    fid.close()
                orivec = arrvec-calctt
                origin_std[iy,ix,iz] = orivec.std()
                origin_mean[iy,ix,iz] = orivec.mean()
    #Now find the 3D index of the best point so far
    min_ind = origin_std.argmin()
    (minx,miny,minz) = search_inds.get_3D(min_ind)
    minx = qx[minx]; miny = qy[miny]; minz = qz[minz];
    return origin_mean,origin_std


class Locator:
    """
    An object class to provide functionality to locate Earthquakes.
    Location parameter configuration is stored in this object class.
    """
    def __init__(self, cfg_dict):
        """
        Initialize locator object with a dictionary of paramaters
        parsed from .cfg file by anfseistools.core.parse_cfg().

        Arguments:
        cfg_dict - Dictionary returned by anfseistools.core.parse_cfg()
        """
        for key in cfg_dict:
            setattr(self, key, cfg_dict[key])

    def locate_eq(self, ev):
        """
        NEEDS TO BE UPDATED
        Locate an earthquake based on the arrivals in ev, traveltime
        files which are already saved.
        """
        loc_params = self.location_parameters
        prop_params = self.propagation_grid
        earth_rad = self.misc['earth_radius']
        #earth_rad = 6371.0


        #Get Propagation grid paramters
        nlat = prop_params['nlat']
        nlon = prop_params['nlon']
        nz = prop_params['nr']
        li  =  LinearIndex(nlon, nlat, nz)
        olon = prop_params['minlon']
        olat = prop_params['minlat']
        oz = prop_params['minz']
        dlon = prop_params['dlon']
        dlat = prop_params['dlat']
        dz = prop_params['dr']

        #Build vectors of geographic coordinates
        qlon = linspace(olon,dlon * nlon + olon,nlon,False)
        qlat = linspace(olat,dlat * nlat + olat,nlat,False)
        qdep = earth_rad - linspace(earth_rad+oz-(nz-1)*dz,earth_rad+oz,nz)
        #print nlat, nlon, nz
        #print len(qlat), len(qlon), len(qdep)
        delta_x = qlon[1] - qlon[0]
        delta_y = qlat[1] - qlat[0]
        delta_z = qdep[1] - qdep[0]

        #Grid search for best location
        start_time = time.time()
        absvec = []
        arrvec = []
        arrsta = []        #a list of station names
        arrpha = []        #List of phases
        for arrival in ev.arrivals:
            if not os.path.isfile('%s%s.traveltime'
                    % (self.misc['tt_map_dir'], arrival.sta)):
                continue
            if arrival.phase is 'P':
                absvec.append(arrival.time)
                arrsta.append(arrival.sta)
                arrpha.append(arrival.phase)
        #if len(arrvec)<1:#Can't locate without at least one
        if len(absvec)<5:#Can't locate without at least one
            print len(ev.arrivals)
            return None
        absvec=asarray(absvec)

        #Search coarsely
        #dstep should go in parameter file.
        dstep = int(loc_params['dstep2'])
        dx, dy, dz = nlon / dstep, nlat / dstep, nz / dstep
        qx, qy, qz = range(1, nlon, dx), range(1, nlat, dy), range(1, nz, dz);
        minx, miny, minz, orgmin,ha = self.grid_search_abs(arrsta, qx, qy,qz, absvec, li)

        #Finer search
        buffx=qx[1]-qx[0]
        buffy=qy[1]-qy[0]
        buffz=buffx/3 #Depth is much trickier to nail down, so this will be done later more strictly
        qx = range(minx - buffx, minx + buffx)
        qy = range(miny - buffy, miny + buffy)
        qz = range(minz - buffz, minz + buffz);
        qx = self.fix_boundary_search(qx, li.nx)
        qy = self.fix_boundary_search(qy, li.ny)
        qz = self.fix_boundary_search(qz, li.nz)
        minx, miny, minz, otime,ha = self.grid_search_abs(arrsta, qx, qy,qz, absvec, li)

        #Get depth
        qx,qy=[minx],[miny]
        qz=range(nz)
        minx, miny, minz, otime,ha = self.grid_search_abs(arrsta, qx, qy,qz, absvec, li)

        #Best-fit grid point
        glon, glat, gz = qlon[minx], qlat[miny], qdep[minz]

        #Get subgrid location
        for i in range(10):#This is really a while loop, but like this in case it is degenerate
            #loc=core_tools.Locator(cfg_dict)
            c,resid,tt_updated,sigma,resid_std=self.get_subgrid_loc(minx,miny,minz,absvec,arrsta,li)
            loc_change=c*[delta_x, delta_y, delta_z]
            #Find the best-fit source location in geographic coordinates
            newloc=[newlon,newlat,newz]=asarray([glon,glat,gz])+loc_change
            ix=nonzero(qlon==find_nearest(qlon,newlon))[0][0]
            iy=nonzero(qlat==find_nearest(qlat,newlat))[0][0]
            iz=nonzero(qdep==find_nearest(qdep,newz))[0][0]
            if minx==ix and miny==iy and minz==iz:
                break
            minx,miny,minz=ix,iy,iz
        #if new location is outside the boundary of the velocity model
        #it can't be trusted, return None
        if newloc[0] < min(qlon) or newloc[0] > max(qlon) or \
                newloc[1] < min(qlat) or newloc[1] > max(qlat) or \
                newloc[2] < min(qdep) or newloc[2] > max(qdep):
            return None
        #Add calculated travel times to Origin
#        ic=0 #Just a counter
#        for arrival in ev.arrivals: #Loop over all input phases
#            if (arrsta[ic] is arrival.sta) and (arrpha[ic] is arrival.phase):
#                arrival.tt_calc=tt_updated[ic] #Only update the ones used in the relocation
#            ic+=1
        for ic in range(len(arrpha)):
            for arrival in ev.arrivals:
                if (arrsta[ic] == arrival.sta) and (arrpha[ic] == arrival.phase):
                    arrival.tt_calc = tt_updated[ic]
                    break

        elapsed_time = time.time() - start_time
        return Origin(newlat,           #Add sigma and residual standard deviation (rsid_std) AAA
                      newlon,
                      newz,
                      otime,
                      'PyLocEQ',
                      arrivals=ev.arrivals,
                      evid=ev.evid,
                      nass=len(ev.arrivals),
                      ndef=len(absvec))

    def grid_search_traveltimes_origin(self, arrsta, qx, qy, qz, arrvec, li):
        """
        NEEDS TO BE UPDATED
        Find the minimum value of the origin time standard deviation following Ben-Zion et al., 1992 (JGR)
           sta          list of station names; strings
           qx,qy,qz    vectors of indices to search through
           arrvec      vector of absolute arrivals in the same order as sta
           li          LinearIndex class for the entire traveltime grid
        """
        from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
        origin_std = array([])
        origin_mean = array([])
        #origin_std = empty([ len(qy),len(qx),len(qz)] )+1000 #Give large starting values
        #origin_mean = empty([ len(qy),len(qx),len(qz)] )
        search_inds = LinearIndex(len(qx),len(qy),len(qz))
        for ix in qx: #range(len(qx)):  #Loop over the three vectors, searching every point
            #print 'On ix: ', ix,'/',len(qx),'\n'
            for iy in qy: #range(len(qy)):
                for iz in qz: #range(len(qz)):
                    calctt = array([]); #initialize the calculated tt vector
                    #ind = li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                    ind = li.get_1D(ix,iy,iz)
                    #Make a traveltime vector from calculated times
                    calctt = read_tt_vector(arrsta,
                                            ind,
                                            self.misc['tt_map_dir'])
                    orivec = arrvec - calctt #Take the difference
                    origin_mean = append( origin_mean,orivec.mean()) #find the mean origin time
                    if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                        origin_std = append(origin_std,10000) #A large dummy value
                        continue
                    origin_std = append(origin_std,orivec.std)
                    #origin_std[iy,ix,iz] = orivec.std()
                    #origin_mean[iy,ix,iz] = orivec.mean()
        #Now find the 3D index of the best point so far
        #new_origin = origin_mean.flatten()[oristd.argmin()] #The origin time with lowest std
        min_ind = origin_std.argmin()
        (minx,miny,minz) = search_inds.get_3D(min_ind)
        minx = qx[minx]; miny = qy[miny]; minz = qz[minz];
        return minx,miny,minz,origin_mean[min_ind]

    def grid_search_traveltimes_rms(self, arrsta, qx, qy, qz, arrvec, li):
        """
        NEEDS TO BE UPDATED
        Find the minimum value of some criterion by performing a grid search
        We aren't necessarily searching the whole grid; we may be skipping values
           sta         list of station names; strings
           qx,qy,qz    vectors of indices to search through
           arrvec      vector of arrivals in the same order as sta
           li          LinearIndex class for the entire traveltime grid
        """
        from numpy import array #There is no reason this should be here, but the next line generated error messages if it wasn't. I'm confused
        rms = array([])
        search_inds = LinearIndex(len(qx),len(qy),len(qz))
        for ix in qx:  #Loop over the three vectors, searching every point
            for iy in qy:
                for iz in qz:
                    calctt = array([]); #initialize the calculated tt vector
                    ind = li.get_1D(ix,iy,iz) #Find the vector index
                    for sta in arrsta: #Build vector of calculated ttimes
                        #print 'bin.'+sta+'.traveltime'
                        fid = open('bin.'+sta+'.traveltime')
                        calctt = append(calctt, read_binary_float(fid,ind) )
                        #print read_binary_float(fid,ind)
                        fid.close()
                    rms = append(rms, sqrt(mean( (arrvec - calctt)**2 )) ) #Root-mean-square, yo
        #Now find the 3D index of the best point so far
        min_ind = rms.argmin()
        (minx,miny,minz) = search_inds.get_3D(min_ind)
        minx = qx[minx]; miny = qy[miny]; minz = qz[minz];
        return minx,miny,minz

    def grid_search_abs(self, arrsta, qx, qy, qz, arrvec, li):
        """
        NEEDS TO BE UPDATED
        Find the minimum of the absolute value of the
                calculated origin time following Ben-Zion et al., 1992 (JGR)
           sta          list of station names; strings
           qx,qy,qz    vectors of indices to search through
           arrvec      vector of absolute arrivals in the same order as sta
           li          Linear_index class for the entire traveltime grid
        """
        from numpy import array,indices #AAA MAY NEED TO DELETE
        best_misfit=100000.0
        search_inds=LinearIndex(len(qx),len(qy),len(qz))
        for ix in range(len(qx)):  #Loop over the three vectors, searching every point
            for iy in range(len(qy)):
                for iz in range(len(qz)):
                    calctt=array([]); #initialize the calculated tt vector
                    ind=li.get_1D(qx[ix],qy[iy],qz[iz]) #Find the vector index
                    calctt=read_tt_vector(arrsta,
                                          ind,
                                          self.misc['tt_map_dir']) #Make traveltime vector from calculated times
                    if min(calctt)<0: #If the traveltime <0, this gridpoint is null
                        continue
                    #Compute misfit measurements
                    #mini=calctt.argmin() #We will compute the origin time from the smallest traveltime Following Pavlis et al., 2004
                    orivec=arrvec-calctt #Take the difference
                    #otime=orivec[mini] #Origin time 
                    otime=orivec.mean()
                    res=orivec-otime
                    #weight=1/abs(res)/(1/abs(res)).sum()  #We can weight by e.g., residual here
                    #misfit=(abs(res)*weight).sum() #Weighted absolute value
                    misfit=(abs(res)).sum() #Absolute value
                    if misfit<best_misfit:
                        best_misfit=misfit
                        minx,miny,minz=qx[ix],qy[iy],qz[iz]
                        best_ot=otime
        return minx,miny,minz,best_ot,best_misfit
        #return minx,miny,minz,origin_mean,misfit

    def get_subgrid_loc(self, ix, iy, iz, arrvec, arrsta, li):
        """
        NEEDS TO BE UPDATED
        """
        #Test least squares on real data
        import numpy as np
        from scipy import linalg

        #Cut off derivative calculations at model boundaries
        endx=ix+1; endy=iy+1; endz=iz+1
        if endx==li.nx: endx=ix
        if endy==li.ny: endy=iy
        if endz==li.nz: endz=iz
        #Get traveltime vectors for the closest point and its neighbors
        ind = li.get_1D(ix,iy,iz)
        tt000 = read_tt_vector(arrsta, ind, self.misc['tt_map_dir'])
        ind = li.get_1D(endx, iy, iz)
        tt100 = read_tt_vector(arrsta, ind, self.misc['tt_map_dir'])
        ind = li.get_1D(ix, endy, iz)
        tt010 = read_tt_vector(arrsta, ind, self.misc['tt_map_dir'])
        ind = li.get_1D(ix, iy, endz)
        tt001 = read_tt_vector(arrsta, ind, self.misc['tt_map_dir'])
        #backwards WILL ALSO NEED EDGEPROOFING !!!!!
        ind=li.get_1D(ix-1,iy,iz)
        btt100= read_tt_vector(arrsta,ind, self.misc['tt_map_dir'])
        ind=li.get_1D(ix,iy-1,iz)
        btt010= read_tt_vector(arrsta,ind, self.misc['tt_map_dir'])
        ind=li.get_1D(ix,iy,iz-1)
        btt001= read_tt_vector(arrsta,ind, self.misc['tt_map_dir'])

        #Calculate forward derivatives
        dt_dx=tt100-tt000
        dt_dy=tt010-tt000
        dt_dz=tt001-tt000
        #backwards
        bdt_dx=tt000-btt100
        bdt_dy=tt000-btt010
        bdt_dz=tt000-btt001
        #Central
        dt_dx=(dt_dx+bdt_dx)/2 #average
        dt_dy=(dt_dy+bdt_dy)/2 #average
        dt_dz=(dt_dz+bdt_dz)/2 #average

        #Build and condition residual
        r=arrvec-tt000
        r=r-r.mean()
        #Create weight based on the sensitivity to depth changes (since it is hardest to constrain)
        w=dt_dz/dt_dz.sum() #Normalize by sum. May be unnecessary
        #Build A matrix, w, and r in wr=Ax  (x is the spatial vector here [x,y,z])
        A=c_[dt_dx,dt_dy,dt_dz]
        wr=r #W*r   Currently not using the weight; doesn't seem necessary here
        c,resid,rank,sigma=linalg.lstsq(A,wr)

        #Compute updated travel times
        tt_updated=tt000+(A*c).sum(axis=1) #c has independent changes for x,y,z, so sum them
        #Compute variance-covariance matrix
        A=c_[dt_dx,dt_dy,dt_dz,dt_dx*0+1] #Add origin time 'derivative'
        sigma = np.dot(A.transpose(),A) #There is probably more to it than this...

        return c, resid,tt_updated,sigma,r.std()

    def fix_boundary_search(self, qx, nx):
        """
        NEEDS TO BE UPDATED
        When performing a grid search on a subgrid, make sure you don't go off the edges
          qx         search vectors, these will be modified then returned
          nx         max index [li.nx]
        """
        for ix in range(len(qx)):
            if qx[ix] < 0:
                qx[ix] = 0
            if qx[ix] >= nx:
                qx[ix] = nx - 1
        newqx = uniq(qx)
        return newqx

def uniq(input):
    """
    NEEDS TO BE UPDATED
    Remove duplicate items from a list. Preserves order.
    """
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

class LinearIndex():
    """
    NEEDS TO BE UPDATED
    Holds a 1D list of 3D indices and a 3D list of 1D indices
    where iz varies fastest, then iy, then ix
    The speed of this can certainly be improved
    """
    def __init__(self, nx, ny, nz):
        from numpy import empty
        self.nx = nx; self.ny = ny; self.nz = nz;
        self.i1D = []
        self.i3D = empty((nx, ny, nz)) #python has weird index conventions
        ic = 0
        for ix in range(nx): #Fuck it, just be explicit
            for iy in range(ny):
                for iz in range(nz):
                    self.i1D.append((ix, iy, iz))
                    self.i3D[ix, iy, iz] = ic
                    ic = ic + 1
        self.i3D = self.i3D.astype(int)

    def get_1D(self, ix, iy, iz):
        return self.i3D[ix, iy, iz]

    def get_3D(self, iv):
        return self.i1D[iv]

class Station:
    """
    A container class for station location data.
    """
    def __init__(self, sta, lat, lon, elev):
        """
        Initialize Station object.

        Arguments:
        sta - Station code.
        lat - Station latitude.
        lon - Station longitude.
        elev - Statio elevation.
        """
        self.sta = sta
        self.lat = lat
        self.lon = lon
        self.elev = elev

    def __str__(self):
        """
        Return string representation of self object.
        """
        ret = 'Station Object\n--------------\n'
        ret += 'sta:\t\t%s\n' % self.sta
        ret += 'lat:\t\t%s\n' % self.lat
        ret += 'lon:\t\t%s\n' % self.lon
        ret += 'elev:\t\t%s\n' % self.elev
        return ret

class Event():
    """
    A container class for Earthquake event data. Mirrors the Event
    table of the CSS3.0 databse schema.
    """
    #def __init__(self, time, lat, lon, depth, mag, magtype=None, evid=None):
    def __init__(self,
                 prefor=None,
                 evid=None,
                 evname=None,
                 auth=None,
                 commid=None,
                 lddate=None,
                 origins=None):
        """
        Initialize Event object.

        Arguments:
        prefor - Preferred origin ID.

        Keyword Arguments:
        evid - Event ID.
        evname - Event author.
        commid - Comment ID.
        lddate - Load date.
        origins - List of anfseistools.core.Origin objects.

        Example:
        In [1]: from anfseistools.core import Event, Origin

        In [2]: origin = Origin(33.4157,
                                -116.8622,
                                4.8910,
                                1275439331.718,
                                orid=287456,
                                nass=47,
                                ndef=47,
                                auth='ANF:vernon',
                                evid=202856,
                                algorithm='locsat:iasp91')

        In [3]: event = Event(prefor=287456,
                              evid=202856,
                              auth='ANF:vernon',
                              origins=[origin])

        In [4]: print event
        Event Object
        ------------
        evid:       202856
        evname:     None
        prefor:     287456
        auth:       ANF:vernon
        commid:     None
        lddate:     None
        origins:
                    Origin Object
                    -------------
                    lat:        33.4157
                    lon:        -116.8622
                    depth:      4.891
                    time:       1275439331.72
                    orid:       287456
                    evid:       202856
                    auth:       ANF:vernon
                    jdate:      None
                    nass:       47
                    ndef:       47
                    ndp:        None
                    grn:        None
                    srn:        None
                    etype:      None
                    review:     None
                    depdp:      None
                    dtype:      None
                    mb:     None
                    mbid:       None
                    ms:     None
                    msid:       None
                    None
                    mlid:       None
                    algorithm:      locsat:iasp91
                    commid:     None
                    lddate:     None
                    arrivals:
        """
        import time as pytime
        self.evid = evid
        self.evname = evname
        self.auth = auth
        self.commid = commid
        self.lddate = lddate
        self.preferred_origin = None
        if origins == None: self.origins = []
        else: self.origins = origins
        self.set_prefor(prefor)

    def __str__(self):
        """
        Return the string representation of anfseistools.core.Event
        object.
        """
        ret = 'Event Object\n------------\n'
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'evname:\t\t%s\n' % self.evname
        ret += 'prefor:\t\t%s\n' % self.prefor
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'origins:\n'
        if len(self.origins) == 0:
            ret += '\t\tNone\n'
        else:
            for i in range(len(self.origins)):
                for line in  ('%s' % self.origins[i]).split('\n'):
                    ret += '\t\t%s\n' % line
        return ret

    def set_prefor(self, prefor):
        """
        Set self.prefor equal to new origin ID and set
        self.preferred_origin to point to the anfseistools.core.Origin
        object referred to by that origin ID.

        Arguments:
        prefor - The origin ID (orid) of the preferred solution.

        Example:
        In [1]: from anfseistools.core import Event, Origin

        In [2]: origin1 = Origin(43.7000,
                                 -79.4000,
                                 5.0,
                                 1398883850.648,
                                 'White',
                                 orid=1234,
                                 evid=1001)

        In [3]: origin2 = Origin(43.7050,
                                 -79.3981,
                                 7.3,
                                 1398883851.346,
                                 'White',
                                 orid=1235,
                                 evid=1001)

        In [4]: event = Event(prefor=1234,
                              evid=1001,
                              auth='White',
                              origins=[origin1, origin2])

        In [5]: print event.preferred_origin
        Origin Object
        -------------
        lat:        43.7
        lon:        -79.4
        depth:      5.0
        time:       1398883850.65
        orid:       1234
        evid:       1001
        auth:       White
        jdate:      None
        nass:       None
        ndef:       None
        ndp:        None
        grn:        None
        srn:        None
        etype:      None
        review:     None
        depdp:      None
        dtype:      None
        mb:     None
        mbid:       None
        ms:     None
        msid:       None
        ml:     None
        mlid:       None
        algorithm:      None
        commid:     None
        lddate:     None
        arrivals:


        In [6]: event.set_prefor(1235)
        Out[6]: 0

        In [7]: print event.preferred_origin
        Origin Object
        -------------
        lat:        43.705
        lon:        -79.3981
        depth:      7.3
        time:       1398883851.35
        orid:       1235
        evid:       1001
        auth:       White
        jdate:      None
        nass:       None
        ndef:       None
        ndp:        None
        grn:        None
        srn:        None
        etype:      None
        review:     None
        depdp:      None
        dtype:      None
        mb:     None
        mbid:       None
        ms:     None
        msid:       None
        ml:     None
        mlid:       None
        algorithm:      None
        commid:     None
        lddate:     None
        arrivals:
        """
        self.prefor = prefor
        for i in range(len(self.origins)):
            if self.origins[i].orid == prefor:
                self.preferred_origin = self.origins[i]
                return 0
        if len(self.origins) == 0:
            return -1
        else:
            self.preferred_origin = self.origins[0]
            return 1

    def add_origin(self,
                   lat,
                   lon,
                   depth,
                   time,
                   auth,
                   arrivals=[],
                   orid=None,
                   evid=None,
                   jdate=None,
                   nass=None,
                   ndef=None,
                   ndp=None,
                   grn=None,
                   srn=None,
                   etype=None,
                   review=None,
                   depdp=None,
                   dtype=None,
                   mb=None,
                   mbid=None,
                   ms=None,
                   msid=None,
                   ml=None,
                   mlid=None,
                   algorithm=None,
                   commid=None,
                   lddate=None):
        """
        Add an anfseistools.core.Origin object to the list of origins
        associated with this event.

        Arguments:
        lat - Latitude of Earthquake hypocenter.
        lon - Longitude of Earthquake hypocenter.
        depth - Depth of Earthquake hypocenter.
        time - Epoch time of Earthquake rupture.

        Keyword Arguments:
        These need to be FULLY described here. Procastinating on this,
        see below.  These fields are optional and exist for posterity
        and to mirror the Origin table of the CSS3.0 schema in whole.
        Refer to CSS3.0 schema for details
        (https://anf.ucsd.edu/pdf/css30.pdf).

        Example:
        In [1]: from anfseistools.core import Event

        In [2]: event = Event(prefor=287456, evid=202856, auth='ANF:vernon')

        In [3]: print event
        Event Object
        ------------
        evid:       202856
        evname:     None
        prefor:     287456
        auth:       ANF:vernon
        commid:     None
        lddate:     None
        origins:
                    None


        In [4]: event.add_origin(33.4157,
                                 -116.8622,
                                 4.8910,
                                 1275439331.718,
                                 orid=287456,
                                 nass=47,
                                 ndef=47,
                                 auth='ANF:vernon',
                                 evid=202856,
                                 algorithm='locsat:iasp91')

        In [5]: print event
        Event Object
        ------------
        evid:       202856
        evname:     None
        prefor:     287456
        auth:       ANF:vernon
        commid:     None
        lddate:     None
        origins:
                    Origin Object
                    -------------
                    lat:        33.4157
                    lon:        -116.8622
                    depth:      4.891
                    time:       1275439331.72
                    orid:       287456
                    evid:       202856
                    auth:       ANF:vernon
                    jdate:      None
                    nass:       47
                    ndef:       47
                    ndp:        None
                    grn:        None
                    srn:        None
                    etype:      None
                    review:     None
                    depdp:      None
                    dtype:      None
                    mb:     None
                    mbid:       None
                    ms:     None
                    msid:       None
                    ml:     None
                    mlid:       None
                    algorithm:      locsat:iasp91
                    commid:     None
                    lddate:     None
                    arrivals:
        """
        self.origins += [Origin(lat,
                                lon,
                                depth,
                                time,
                                auth,
                                orid=orid,
                                evid=evid,
                                arrivals=arrivals,
                                jdate=jdate,
                                nass=nass,
                                ndef=ndef,
                                ndp=ndp,
                                grn=grn,
                                srn=srn,
                                etype=etype,
                                review=review,
                                depdp=depdp,
                                dtype=dtype,
                                mb=mb,
                                mbid=mbid,
                                ms=ms,
                                msid=msid,
                                ml=ml,
                                mlid=mlid,
                                algorithm=algorithm,
                                commid=commid,
                                lddate=lddate)]
class Origin():
    """
    A container class for Earthquake event data. Mirrors the Origin
    table of the CSS3.0 databse schema.
    """
    def __init__(self,
                 lat,
                 lon,
                 depth,
                 time,
                 auth,
                 arrivals=[],
                 orid=None,
                 evid=None,
                 jdate=None,
                 nass=None,
                 ndef=None,
                 ndp=None,
                 grn=None,
                 srn=None,
                 etype=None,
                 review=None,
                 depdp=None,
                 dtype=None,
                 mb=None,
                 mbid=None,
                 ms=None,
                 msid=None,
                 ml=None,
                 mlid=None,
                 algorithm=None,
                 commid=None,
                 lddate=None):
        """
        Initialize Origin object.

        Arguments:
        lat - Latitude of Earthquake hypocenter.
        lon - Longitude of Earthquake hypocenter.
        depth - Depth of Earthquake hypocenter.
        time - Epoch time of Earthquake rupture.

        Keyword Arguments:
        These need to be FULLY described here. Procastinating on this,
        see below.  These fields are optional and exist for posterity
        and to mirror the Origin table of the CSS3.0 schema in whole.
        Refer to CSS3.0 schema for details
        (https://anf.ucsd.edu/pdf/css30.pdf).

        Example:
        In [1]: from anfseistools.core import Origin, Arrival

        In [2]: arrivals = [Arrival('SND',
                                  1276817657.230,
                                  'P',
                                  chan='HHZ')]

        In [3]: arrivals += [Arrival('FRD',
                                   1276817656.000,
                                   'P',
                                   chan='HHZ')]

        In [4]: origin = Origin(32.7103,
                                -115.9378,
                                3.44,
                                1276817637.470,
                                'White',
                                orid=235993,
                                evid=2010168,
                                nass=6,
                                ndef=64)

        In [5]: print origin
        Origin Object
        -------------
        lat:        32.7103
        lon:        -115.9378
        depth:      3.44
        time:       1276817637.47
        orid:       235993
        evid:       2010168
        auth:       White
        jdate:      None
        nass:       6
        ndef:       64
        ndp:        None
        grn:        None
        srn:        None
        etype:      None
        review:     None
        depdp:      None
        dtype:      None
        mb:     None
        mbid:       None
        ms:     None
        msid:       None
        ml:     None
        mlid:       None
        algorithm:      None
        commid:     None
        lddate:     None
        arrivals:
        """

        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.time = time
        self.orid = orid
        self.evid = evid
        self.auth = auth
        self.arrivals = arrivals
        self.jdate = jdate
        self.nass = nass
        self.ndef = ndef
        self.ndp = ndp
        self.grn = grn
        self.srn = srn
        self.etype = etype
        self.review = review
        self.depdp = depdp
        self.dtype = dtype
        self.mb = mb
        self.mbid = mbid
        self.ms = ms
        self.msid = msid
        self.ml = ml
        self.mlid = mlid
        self.algorithm = algorithm
        self.commid = commid
        self.lddate = lddate

    def __str__(self):
        """
        Returns the string representation of anfseistools.core.Origin
        object.
        """
        ret = 'Origin Object\n-------------\n'
        ret += 'lat:\t\t%s\n' % self.lat
        ret += 'lon:\t\t%s\n' % self.lon
        ret += 'depth:\t\t%s\n' % self.depth
        ret += 'time:\t\t%s\n' % self.time
        ret += 'orid:\t\t%s\n' % self.orid
        ret += 'evid:\t\t%s\n' % self.evid
        ret += 'auth:\t\t%s\n' % self.auth
        ret += 'jdate:\t\t%s\n' % self.jdate
        ret += 'nass:\t\t%s\n' % self.nass
        ret += 'ndef:\t\t%s\n' % self.ndef
        ret += 'ndp:\t\t%s\n' % self.ndp
        ret += 'grn:\t\t%s\n' % self.grn
        ret += 'srn:\t\t%s\n' % self.srn
        ret += 'etype:\t\t%s\n' % self.etype
        ret += 'review:\t\t%s\n' % self.review
        ret += 'depdp:\t\t%s\n' % self.depdp
        ret += 'dtype:\t\t%s\n' % self.dtype
        ret += 'mb:\t\t%s\n' % self.mb
        ret += 'mbid:\t\t%s\n' % self.mbid
        ret += 'ms:\t\t%s\n' % self.ms
        ret += 'msid:\t\t%s\n' % self.msid
        ret += 'ml:\t\t%s\n' % self.ml
        ret += 'mlid:\t\t%s\n' % self.mlid
        ret += 'algorithm:\t\t%s\n' % self.algorithm
        ret += 'commid:\t\t%s\n' % self.commid
        ret += 'lddate:\t\t%s\n' % self.lddate
        ret += 'arrivals:\n'
        for i in range(len(self.arrivals)):
            ret += '%s' % self.arrivals[i]
        return ret

    def update_predarr_times(self, cfg_dict):
        """
        Update the anfseistools.core.Arrival.tt_calc and
        anfseistools.core.Arrival.predarr fields for each Arrival object in
        anfseistools.core.Origin object's arrivals attribute.

        Arguments:
        cfg_dict - Configuration dictionary as returned by
        anfseistools.core.parse_cfg()

        Caveats:
        This functionality is currently only implemented for P-wave
        arrivals.

        Example:
        In [1]: import sys

        In [2]: import os

        In [3]: sys.path.append('%s/data/python' % os.environ['ANTELOPE'])

        In [4]: from antelope.datascope import closing, dbopen

        In [5]: import anfseistools.core as core

        In [6]: import anfseistools.ant as ant

        In [7]: cfg_dict = core.parse_cfg('test_pf_2_cfg.cfg')

        In [8]: locator = core.Locator(cfg_dict)

        In [9]: with closing(dbopen('/Users/mcwhite/staging/dbs/' \
                                    'June2010/June2010', 'r')) as db:
           ...:     tbl_event = db.schema_tables['event']
           ...:     tbl_event = tbl_event.subset('evid == 202856')
           ...:     events = ant.create_event_list(tbl_event)
           ...:

       In [10]: new_origin = locator.locate_eq(events[0].preferred_origin)

       In [11]: for arrival in new_origin.arrivals:
          ....:     print arrival.phase, arrival.predarr
          ....:
        S None
        P None
        S None
        P None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        S None
        S None
        P None
        S None
        S None
        S None
        S None
        S None
        S None
        P None
        S None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        P None
        S None
        S None
        S None
        P None
        S None
        P None

        In [12]: new_origin.update_predarr_times(cfg_dict)
        Out[12]: 0

        In [13]: for arrival in new_origin.arrivals:
            ...:     print arrival.phase, arrival.predarr
            ...:
        S None
        P 1275439337.07
        S None
        P 1275439337.66
        P 1275439338.11
        S None
        P 1275439338.22
        S None
        P 1275439338.29
        S None
        P 1275439339.54
        S None
        S None
        S None
        P 1275439339.66
        S None
        S None
        S None
        S None
        S None
        S None
        P 1275439344.33
        S None
        S None
        P 1275439337.06
        S None
        P 1275439335.57
        S None
        P 1275439335.58
        S None
        P 1275439335.84
        S None
        P 1275439336.15
        S None
        P 1275439336.19
        S None
        P 1275439336.4
        S None
        P 1275439336.6
        S None
        P 1275439336.6
        S None
        S None
        S None
        P 1275439337.19
        S None
        P 1275439333.59
        """
        from numpy import linspace
        #Get Propagation grid paramters
        ttdir = cfg_dict['misc']['tt_map_dir']
        prop_params = cfg_dict['propagation_grid']
        earth_rad = cfg_dict['misc']['earth_radius']
        nlat = prop_params['nlat']
        nlon = prop_params['nlon']
        nz = prop_params['nr']
        li  =  LinearIndex(nlon, nlat, nz)
        olon = prop_params['minlon']
        olat = prop_params['minlat']
        oz = prop_params['minz']
        dlon = prop_params['dlon']
        dlat = prop_params['dlat']
        dz = prop_params['dr']
        #Build vectors of geographic coordinates
        qlon = linspace(olon, dlon * nlon + olon, nlon, False)
        qlat = linspace(olat, dlat * nlat + olat, nlat, False)
        qdep = earth_rad - linspace(earth_rad + oz - (nz - 1) * dz,
                                    earth_rad + oz,
                                    nz)
        endpoints, indices = find_containing_cube(self.lat,
                                                  self.lon,
                                                  self.depth,
                                                  qlat,
                                                  qlon,
                                                  qdep)
        for arrival in self.arrivals:
            if arrival.phase == 'P':
                if not os.path.isfile('%sbin.%s.traveltime'
                        % (ttdir, arrival.sta)):
                    continue
                ttvec = []
                for i in range(len(indices)):
                    index, endpoint = indices[i], endpoints[i]
                    li1D = li.get_1D(index[1], index[0], index[2])
                    ttvec += [read_tt_vector([arrival.sta], li1D, ttdir)[0]]
                dtt_dlat =  0 if endpoints[1][0] == endpoints[0][0] else\
                        (ttvec[1] - ttvec[0]) / \
                        (endpoints[1][0] - endpoints[0][0])
                dtt_dlon = 0 if endpoints[3][1] == endpoints[0][1] else\
                        (ttvec[3] - ttvec[0]) / \
                        (endpoints[3][1] - endpoints[0][1])
                dtt_ddep = 0 if endpoints[4][2] == endpoints[0][2] else\
                        (ttvec[4] - ttvec[0]) / \
                        (endpoints[4][2] - endpoints[0][2])
                delta_lon = self.lon - endpoints[0][1]
                delta_lat = self.lat - endpoints[0][0]
                delta_dep = self.depth - endpoints[0][2]
                tt = ttvec[0] + (dtt_dlon * delta_lon)\
                            + (dtt_dlat * delta_lat)\
                            + (dtt_ddep * delta_dep)
                predarr = self.time + tt
                arrival.tt_calc = tt
                arrival.predarr = predarr
        return 0

class Arrival():
    """
    A container class for phase data.
    """
    def __init__(self,
                 sta,
                 time,
                 phase,
                 chan=None,
                 deltim=None,
                 qual=None,
                 arid=None):
        """
        Initialize anfseistools.core.Arrival object.

        Arguments:
        sta - Station name.
        time - Epoch time of arrival observation.
        phase - Phase type (Eg. P, Pn, Pb, Pg, S, Sn, Sb, Sg)
        chan - Channel observation made on.
        deltim - Standard deviation of observed arrival time.
        qual - Signal onset quality
        (Ie. i: impulsive, e: emergent, w: weak).
        arid - Arrival ID.
        """
        self.sta = sta
        self.time = time
        self.phase = phase
        self.chan = chan
        self.deltim = deltim
        self.qual = qual
        self.arid = arid
        self.tt_calc = None #calculated travel time
        self.predarr = None #predicted arrival time 

    def __str__(self):
        """
        Return string representation for anfseistools.core.Arrival
        object.
        """
        ret = 'Arrival Object\n--------------\n'
        ret += 'sta:\t\t%s\n' % self.sta
        ret += 'time:\t\t%s\n' % self.time
        ret += 'phase:\t\t%s\n' % self.phase
        ret += 'arid:\t\t%s\n' % self.arid
        ret += 'deltim:\t\t%s\n' % self.deltim
        ret += 'qual:\t\t%s\n'  % self.qual
        ret += 'tt_calc:\t\t%s\n' % self.tt_calc
        ret += 'predarr:\t\t%s\n' % self.predarr
        return ret

import sys
import os
import logging
import time
import subprocess
import shutil
from anfseistools.core import Station, num
from math import pi

tt_calculator = 'fm3d'
tt_dir = 'tt_maps_%d' % int(time.time())
get_time = lambda: time.strftime('%m/%d/%Y %H:%M:%S')

def _parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('db', type=str, help='Input database.')
    parser.add_argument('-b', '--binary', action='store_true',
        help='Write binary output file.')
    parser.add_argument('-p', '--pf', type=str, help='Parameter file.')
    parser.add_argument('-s', '--subset', type=str, help='Station subset.')
    return parser.parse_args()

def _parse_pfile(pfile):
    from antelope.stock import pfin
    return eval_pfile(pfin(pfile).pf2dict())

def eval_pfile(pfile):
    """
    Recursively typecast dictionary values of str-type to int-type
    or float-type values if appropriate.

    The method antelope.stock.ParameterFile.pf2dict returns a
    dictionary with all str-type values. This method is intended
    primarily to take such a dictionary and typecast integers to
    int-type values and floating points to float-type values.

    Arguments:
    pfile - Dictionary to be typecast.

    Returns:
    pfile - Typecasted dictionary.
    """
    for key in pfile:
        if isinstance(pfile[key], dict):
            eval_pfile(pfile[key])
        else:
            if key in locals():
                continue
            try:
                pfile[key] = eval(pfile[key])
            except (NameError, SyntaxError):
                pass

    return pfile

def _check_pfile(pfile):
    exit_flag = False
    prop_grid = pfile['propagation_grid']
    vgrids_infile = open(pfile['vgrids.in'], 'r')
    vgrids_infile.readline()
    vgrids_nr, vgrids_nlat, vgrids_nlon = vgrids_infile.readline().split()
    vgrids_dr, vgrids_dlat, vgrids_dlon = vgrids_infile.readline().split()
    vgrids_minr, vgrids_minlat, vgrids_minlon = vgrids_infile.readline().split()
    vgrids_infile.close()
    vgrids_nr = int(vgrids_nr)
    vgrids_nlat = int(vgrids_nlat)
    vgrids_nlon = int(vgrids_nlon)
    vgrids_dr = float(vgrids_dr)
    vgrids_dlat = float(vgrids_dlat) * 180.0 / pi
    vgrids_dlon = float(vgrids_dlon) * 180.0 / pi
    vgrids_minr = float(vgrids_minr)
    vgrids_minlat = float(vgrids_minlat) * 180.0 / pi
    vgrids_minlon = float(vgrids_minlon) * 180.0 / pi
    vgrids_maxr = vgrids_minr + (vgrids_nr - 1) * vgrids_dr
    vgrids_maxlat = vgrids_minlat + (vgrids_nlat - 1) * vgrids_dlat
    vgrids_maxlon = vgrids_minlon + (vgrids_nlon - 1) * vgrids_dlon
    vgrids_violation_msg = "Propagation grid must be entirely contained by the"\
            " second layer of the velocity grid. Edit the parameter file to "\
            "satisfy the below requirement(s):"
    print_vgrids_violation_msg = True
    if prop_grid['minr'] < (vgrids_minr + vgrids_dr):
        exit_flag = True
        if print_vgrids_violation_msg:
            print vgrids_violation_msg
            print_vgrids_violation_msg = False
        print "\t%.4f <= propagation_grid['minr'] <= %.4f" \
                % ((vgrids_minr + vgrids_dr), (vgrids_maxr - vgrids_dr))
    if prop_grid['minr'] + (prop_grid['nr'] - 1) * prop_grid['dr'] > \
                (vgrids_maxr - vgrids_dr):
        exit_flag = True
        if print_vgrids_violation_msg:
            print vgrids_violation_msg
            print_vgrids_violation_msg = False
        print "\t%.4f <= propagation_grid['minr']  + "\
                "(propagation_grid['nr'] - 1) * propagation_grid['dr'] <= "\
                "%.4f" % ((vgrids_minr + vgrids_dr), (vgrids_maxr - vgrids_dr))
    if prop_grid['minlat'] < (vgrids_minlat + vgrids_dlat):
        exit_flag = True
        if print_vgrids_violation_msg:
            print vgrids_violation_msg
            print_vgrids_violation_msg = False
        print "\t%.4f <= propagation_grid['minlat'] <= %.4f" \
                % ((vgrids_minlat + vgrids_dlat), (vgrids_maxlat - vgrids_dlat))
    if prop_grid['minlat'] + (prop_grid['nlat'] - 1) * prop_grid['dlat'] > \
                (vgrids_maxlat - vgrids_dlat):
        exit_flag = True
        if print_vgrids_violation_msg:
            print vgrids_violation_msg
            print_vgrids_violation_msg = False
        print "\t%.4f <= propagation_grid['minlat']  + "\
                "(propagation_grid['nlat'] - 1) * propagation_grid['dlat'] <= "\
                "%.4f" % ((vgrids_minlat + vgrids_dlat),
                          (vgrids_maxlat - vgrids_dlat))
    if prop_grid['minlon'] < (vgrids_minlon + vgrids_dlon):
        exit_flag = True
        if print_vgrids_violation_msg:
            print vgrids_violation_msg
            print_vgrids_violation_msg = False
        print "\t%.4f <= propagation_grid['minlon'] <= %.4f" \
                % ((vgrids_minlon + vgrids_dlon), (vgrids_maxlon - vgrids_dlon))
    if prop_grid['minlon'] + (prop_grid['nlon'] - 1) * prop_grid['dlon'] > \
                (vgrids_maxlon - vgrids_dlon):
        exit_flag = True
        if print_vgrids_violation_msg:
            print vgrids_violation_msg
            print_vgrids_violation_msg = False
        print "\t%.4f <= propagation_grid['minlon']  + "\
                "(propagation_grid['nlon'] - 1) * propagation_grid['dlon'] <= "\
                "%.4f" % ((vgrids_minlon + vgrids_dlon),
                          (vgrids_maxlon - vgrids_dlon))
    interfaces_infile = open(pfile['interfaces.in'], 'r')
    interfaces_infile.readline()
    interfaces_nlat, interfaces_nlon = interfaces_infile.readline().split()
    interfaces_dlat, interfaces_dlon = interfaces_infile.readline().split()
    interfaces_minlat, interfaces_minlon = interfaces_infile.readline().split()
    interfaces_infile.close()
    interfaces_nlat = int(interfaces_nlat)
    interfaces_nlon = int(interfaces_nlon)
    interfaces_dlat = float(interfaces_dlat) * 180.0 / pi
    interfaces_dlon = float(interfaces_dlon) * 180.0 / pi
    interfaces_minlat = float(interfaces_minlat) * 180.0 / pi
    interfaces_minlon = float(interfaces_minlon) * 180.0 / pi
    interfaces_maxlat = interfaces_minlat + (interfaces_nlat - 1) \
            * interfaces_dlat
    interfaces_maxlon = interfaces_minlon + (interfaces_nlon - 1) \
            * interfaces_dlon
    interfaces_violation_msg = "Propagation grid must be entirely contained "\
            "by the second layer of the interfaces grid. Edit the parameter "\
            "file to satisfy the below requirement(s):"
    print_interfaces_violation_msg = True
    if prop_grid['minlat'] < (interfaces_minlat + interfaces_dlat):
        exit_flag = True
        if print_interfaces_violation_msg:
            print interfaces_violation_msg
            print_interfaces_violation_msg = False
        print "\t%.4f <= propagation_grid['minlat'] <= %.4f" \
                % ((interfaces_minlat + interfaces_dlat),
                   (interfaces_maxlat - interfaces_dlat))
    if prop_grid['minlat'] + (prop_grid['nlat'] - 1) * prop_grid['dlat'] > \
                (interfaces_maxlat - interfaces_dlat):
        exit_flag = True
        if print_interfaces_violation_msg:
            print interfaces_violation_msg
            print_interfaces_violation_msg = False
        print "\t%.4f <= propagation_grid['minlat']  + "\
                "(propagation_grid['nlat'] - 1) * propagation_grid['dlat'] <= "\
                "%.4f" % ((interfaces_minlat + interfaces_dlat),
                          (interfaces_maxlat - interfaces_dlat))
    if prop_grid['minlon'] < (interfaces_minlon + interfaces_dlon):
        exit_flag = True
        if print_interfaces_violation_msg:
            print interfaces_violation_msg
            print_interfaces_violation_msg = False
        print "\t%.4f <= propagation_grid['minlon'] <= %.4f" \
                % ((interfaces_minlon + interfaces_dlon),
                   (interfaces_maxlon - interfaces_dlon))
    if prop_grid['minlon'] + (prop_grid['nlon'] - 1) * prop_grid['dlon'] > \
                (interfaces_maxlon - interfaces_dlon):
        exit_flag = True
        if print_interfaces_violation_msg:
            print interfaces_violation_msg
            print_interfaces_violation_msg = False
        print "\t%.4f <= propagation_grid['minlon']  + "\
                "(propagation_grid['nlon'] - 1) * propagation_grid['dlon'] <= "\
                "%.4f" % ((interfaces_minlon + interfaces_dlon),
                          (interfaces_maxlon - interfaces_dlon))
    if exit_flag:
        sys.exit(-1)

def _create_station_list(args, pfile):
    from antelope.datascope import closing, dbopen
    dbpath = args.db
    subset = args.subset
    prop_grid = pfile['propagation_grid']
    minlon = float(prop_grid['minlon'])
    dlon = float(prop_grid['dlon'])
    nlon = int(prop_grid['nlon'])
    minlat = float(prop_grid['minlat'])
    dlat = float(prop_grid['dlat'])
    nlat = int(prop_grid['nlat'])
    maxlon = minlon + dlon * (nlon - 1)
    maxlat = minlat + dlat * (nlat - 1)
    station_list = []
    with closing(dbopen(dbpath, 'r')) as db:
        tbl_site = db.schema_tables['site']
        if subset: tbl_site = tbl_site.subset(subset)
        subset = 'lat > %f && lat < %f && lon > %f && lon < %f' \
                % (minlat, maxlat, minlon, maxlon)
        tbl_site = tbl_site.subset(subset)
        for record in tbl_site.iter_record():
            sta, lat, lon, elev = record.getv('sta',
                                              'lat',
                                              'lon',
                                              'elev')
            station_list += [Station(sta, lat, lon, elev)]
    return station_list

def _write_propgrid(pfile):
    prop_grid = pfile['propagation_grid']
    if os.path.isfile('propgrid.in'):
        os.remove('propgrid.in')
    outfile = open('propgrid.in', 'w')
    outfile.write('%8s\t%8s\t%8s\t\t\t%s\n'
            % (prop_grid['nr'],
               prop_grid['nlat'],
               prop_grid['nlon'],
               '# of propagation grid points in r, lat and lon'))
    outfile.write('%8s\t%8s\t%8s\t\t%s\n'
            % (prop_grid['dr'],
               prop_grid['dlat'],
               prop_grid['dlon'],
               'grid intervals in r (km) lat,long (deg)'))
    outfile.write('%8s\t%8s\t%8s\t\t%s\n'
            % ((prop_grid['minr'] + (prop_grid['nr'] - 1)  * prop_grid['dr'] \
                    - pfile['earth_radius']),
               prop_grid['minlat'],
               prop_grid['minlon'],
               'origin of the grid height (km),lat,long (deg)'))
    outfile.write('%8s\t%8s\t\t\t\t\t%s\n'
            % (prop_grid['refinement_factor'],
               prop_grid['ncells'],
               'refinement factor and # of propgrid cells in refined source '\
                    'grid'))
    outfile.close()

def _configure_logger():
    logger = logging.getLogger(sys.argv[0])
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(fmt='%(asctime)s::%(levelname)s::%(message)s',
        datefmt='%Y%j %H:%M:%S')
    fh = logging.FileHandler('%s.log' % sys.argv[0])
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

def gen_sta_tt_maps(stalist, if_write_binary=True):
#Generate a travel time map for each station in the station list
    print '%s::Starting travel time calculation.' \
            % get_time()
    for sta in stalist:
        print '%s::Generating travel times for station %s.' \
                % (get_time(), sta.sta)
        #Elevation can be set to a large negative number to glue the source to the surface
        #_write_sources_file(sta.elev * -1, sta.lat, sta.lon) #Elevation is in km and negative
        _write_sources_file(0.0, sta.lat, sta.lon) #!!!! SET SOURCE TO 0 DEPTH
        tt_calc_subprocess = subprocess.Popen(tt_calculator,
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.STDOUT)
        for line in tt_calc_subprocess.stdout:
            print '%s::%s' % (get_time(), line.rstrip())
        #Create output file name
        outfnam='%s.traveltime' % sta.sta
        try:
            print '%s::Moving arrtimes.dat to %s' % (get_time(), outfnam)
            shutil.move('arrtimes.dat', outfnam)
        except IOError:
            print '%s::Could not move arrtimes.dat to %s' % (get_time(), outfnam)
            sys.exit(-1)
        if if_write_binary:
            _tt_ascii_2_binary(outfnam)
    elapsed_time=time.time()-start_time
    print '%s::Finished travel time calculation.' % get_time()

def _generate_tt_maps(db, write_binary=True):
    logger = logging.getLogger(sys.argv[0])
    logger.debug('Begin travel-time map generation.')
    with closing(dbopen(db, 'r')) as db:
        tbl_site = db.schema_tables['site']
        for record in tbl_site.iter_record():
            sta, lat, lon, elev = record.getv('sta', 'lat', 'lon', 'elev')
            logger.debug('Begin travel-time map generation for station %s'
                % sta)
            _write_sources_file(0.10, lat, lon)
            os.system(tt_calculator)
            logger.debug('End travel-time map generation for station %s'
                % sta)
#            outfile = '%s.traveltime' % sta
#            os.system('mv arrtimes.dat %s' % outfile)
#            if write_binary:
#                logger.debug('Begin writing binary file for station % s' % sta)
#                _tt_ascii_2_binary(outfile)
#                logger.debug('End writing binary file for station % s' % sta)

def _write_sources_file(depth, lat, lon):
    #Write  the sources.in file, which has this form:
    # 1                                number of sources
    # 0                                source is local/teleseismic (0/1)
    # 5.00  33.0   -116.0      position depth(km),lat(deg),long(deg) 
    # 1                                number of paths from this source
    # 1                                number of sections on the path
    # 0 1           define the path sections 
    # 1            define the velocity type along the path
    numsrc=1
    teleflag=0
    numpaths=1
    numsections=1
    veltype=1
    outfile = open('sources.in','w')
    outfile.write(" %i\n %i\n %.4f %.4f %.4f\n %i\n %i\n %i %i\n %i\n"
        % (numsrc, teleflag, depth, lat, lon, numpaths, numsections, 0, 2,
        veltype))
    outfile.close()
    #Path sections are described in more detail in the FMM README
    # For the first arrival from a source propagating at the surface [ 0 2 ]
    #  is the path section

def _tt_ascii_2_binary(fnam):
#Convert an ascii traveltime file to binary format
# Also puts the header in a separate ascii file
# just put "bin" and "hdr" in front of the filename
    from array import array
    bin_path = 'bin.%s' % fnam
    hdr_path = 'hdr.%s' % fnam
    print '%s::Reading arrival times from %s' % (get_time(), fnam)
    #Open and read header
    fid = open('%s' % fnam, 'r')
    tmp=fid.readline().strip().split()
    nz=num(tmp[0]);nlat=num(tmp[1]);nlon=num(tmp[2]) #Number of grid points
    tmp=fid.readline().strip().split()
    dz=num(tmp[0]);dlat=num(tmp[1]);dlon=num(tmp[2]) #Spacing of grid points
    tmp=fid.readline().strip().split()
    oz=num(tmp[0]);olat=num(tmp[1]);olon=num(tmp[2]) #Origin of grid points
    tmp=fid.readline().strip().split()
    narr=num(tmp[0])                         #number of sets of arrival times
    tmp=fid.readline().strip().split()
    null1=num(tmp[0]);null2=num(tmp[1]);null3=num(tmp[2])                         #source and path for arrival time ???
    #Now read the traveltimes into a list
    lines=fid.readlines()
    data=[]
    for line in lines:
        data.append( float(line) )
    fid.close()
    # Now output in binary format
    fid=open(bin_path,'w')
    out_array=array('d',data)#It is apparently faster to turn this into an array
    out_array.tofile(fid)
    fid.close()

    #Now Write the header
    fid = open(hdr_path,'w')
    fid.write(" %i %i %i\n %.4f %.6f %.6f\n %.5f %.5f %.5f\n %i\n %i %i %i" % (nz,nlat,nlon, dz,dlat,dlon, oz,olat,olon, narr, null1,null2,null3))
    fid.close()
    print '%s::Finished writing arrival times to %s and %s' \
            % (get_time(), bin_path, hdr_path)
    #return data #Don't forget to remove this!!

def gen_tt_map_malcolm(db):
    from antelope.datascope import closing, dbopen
    station_list = StationList(db, is_db=True)
    gen_sta_tt_maps(station_list)

def read_arrtimes(fnam='arrtimes.dat',ifplot=0): #default name
#Read arrival times from file fnam
    print '%s::Reading arrival times from %s' % (get_time(), fnam)
    #Open and read header
    fid = open(fnam,'r')
    tmp=fid.readline().strip().split()
    nz=num(tmp[0]);nlat=num(tmp[1]);nlon=num(tmp[2]) #Number of grid points
    tmp=fid.readline().strip().split()
    dz=num(tmp[0]);dlat=num(tmp[1]);dlon=num(tmp[2]) #Spacing of grid points
    tmp=fid.readline().strip().split()
    oz=num(tmp[0]);olat=num(tmp[1]);olon=num(tmp[2]) #Origin of grid points
    tmp=fid.readline().strip().split()
    poo=num(tmp[0])                         #number of sets of arrival times
    tmp=fid.readline().strip().split()
    poo=num(tmp[0])                         #source and path for arrival time ???

    #Now read the traveltimes and sort into a matrix
    data=empty((nlat,nlon,nz))
    for ix in range(nlon):
        for iy in range(nlat):
            for iz in range(nz):
                tmp=float(fid.readline().strip())
                data[iy,ix,iz]=tmp
    fid.close()

    #Create vectors of geographic coordinates
    elon=olon+(dlon*nlon);
    elat=olat+(dlat*nlat);
    lonvec=linspace(olon,elon,nlon);
    latvec=linspace(olat,elat,nlat);

    #Now plot
    if ifplot: #adhoc
        flon,flat=load_faults()
        plt.figure(figsize=(6,6))
        plt.plot(flon,flat,'k')
        plt.axis([-118,-115,32,35])
        plt.pcolor(lonvec,latvec,data[:,:,6])
    return data

class TomoDD_MOD():
    #Contains the velocity and metadata from a TomoDD MOD file which looks like
    #bld NX NY NZ
    # X1 X2... XN
    # Y1 Y2... YN
    # Z1 Z2... ZN
    #V111 V211 V311 ... VN11
    #V121 V221 V321 ... VN21
    #...
    #V1N1 V2N1 V3N1 ... VNN1
    #V112 V212 V312 ... VN12
    #...
    #V11N V21N V31N ... VN1N
    #...
    # THEN, the entire V matrix is repeated for Vp/Vs
    # so, total file length is NY*NZ*2 + 4
    #
    # NOTE: Left-handed coordinates; depth is positive
    def __init__(self):
        pass
    def read(self,fnam='MOD'):
        #Read a TomoDD MOD velocity file. It is an ascii file that looks like:
        #Open and read header
        self.fnam=fnam
        fid = open(fnam,'r')
        tmp=fid.readline().strip().split()
        self.bld=num(tmp[0])
        self.nx=num(tmp[1])
        self.ny=num(tmp[2])
        self.nz=num(tmp[3]) #Number of grid points
        tmp=fid.readline().strip().split()
        self.qx=asarray([float(ii) for ii in tmp]) 
        tmp=fid.readline().strip().split()
        self.qy=asarray([float(ii) for ii in tmp])
        tmp=fid.readline().strip().split()
        self.qz=asarray([float(ii) for ii in tmp])
        #Loop to create a numpy matrix
        self.vel=empty((self.nx,self.ny,self.nz)) #Check that this imports properly from numpy
        for iz in range(self.nz):
            for iy in range(self.ny):
                tmp=fid.readline().strip().split()
                for ix in range(self.nx):
                    self.vel[ix,iy,iz]=float(tmp[ix])
        self.velrat=empty((self.nx,self.ny,self.nz))
        for iz in range(self.nz):
            for iy in range(self.ny):
                tmp=fid.readline().strip().split()
                for ix in range(self.nx):
                    self.velrat[ix,iy,iz]=float(tmp[ix])
        fid.close()
    def write(self,outfnam='out.MOD'): #Write to MOD format
        self.outfnam=outfnam
        fid=open(outfnam,'w')
        outs=str(self.bld)+' '+str(self.nx)+' '+str(self.ny)+' '+str(self.nz)+'\n'
        fid.write(outs)
        outs=''.join(str(ii)+' ' for ii in self.qx)+'\n'
        fid.write(outs)
        outs=''.join(str(ii)+' ' for ii in self.qy)+'\n'
        fid.write(outs)
        outs=''.join(str(ii)+' ' for ii in self.qz)+'\n'
        fid.write(outs)
        for iz in range(self.nz):
            for iy in range(self.ny):
                outs=''.join('{1:.{0}f} '.format(3,ii)  for ii in self.vel[:,iy,iz])+'\n'
                fid.write(outs)
        for iz in range(self.nz):
            for iy in range(self.ny):
                outs=''.join('{1:.{0}f} '.format(3,ii)  for ii in self.velrat[:,iy,iz])+'\n'
                fid.write(outs)
        fid.close()

def _create_dummy_receivers_file(pfile):
    outfile = open('receivers.in', 'w')
    outfile.write('1\n')
    prop_grid = pfile['propagation_grid']
    minlat = prop_grid['minlat']
    minlon = prop_grid['minlon']
    maxz = pfile['earth_radius'] - prop_grid['minr']
    maxlat = minlat + (prop_grid['nlat'] - 1) * prop_grid['dlat']
    maxlon = minlon + (prop_grid['nlon'] - 1) * prop_grid['dlon']
    minz = maxz - (prop_grid['nr'] - 1) * prop_grid['dr']
    midlat = (minlat + maxlat) / 2.0
    midlon = (minlon + maxlon) / 2.0
    midz = (minz + maxz) / 2.0
    outfile.write('%f    %f    %f\n1\n1\n1' % (midz, midlat, midlon))
    outfile.close()

class Traveltime_header_file():
#This reads the header file passed as an argument by fnam, saves all of the header info to a class
    def __init__(self,fnam):
       fid = open(fnam,'r')
       tmp=fid.readline().strip().split()
       self.nz=num(tmp[0]);self.nlat=num(tmp[1]);self.nlon=num(tmp[2]) #Number of grid points
       tmp=fid.readline().strip().split()
       self.dz=num(tmp[0]);self.dlat=num(tmp[1]);self.dlon=num(tmp[2]) #Spacing of grid points
       tmp=fid.readline().strip().split()
       self.oz=num(tmp[0]);self.olat=num(tmp[1]);self.olon=num(tmp[2]) #Origin of grid points
       tmp=fid.readline().strip().split()
       self.arrsets=num(tmp[0])                         #number of sets of arrival times
       tmp=fid.readline().strip().split()
       self.srcpath=num(tmp[0])                         #source and path for arrival time ???
       fid.close()

if __name__ == '__main__':
    args = _parse_args()
    try:
        print "%s::Creating working directory - %s." % (get_time(), tt_dir)
        os.mkdir(tt_dir)
    except OSError as err:
        print "%s::Could not create working directory. Make sure you have "\
                "write permission  for %s\n\t%s" % (get_time(), os.getcwd(), err)
        sys.exit(-1)
    try:
        os.chdir(tt_dir)
    except OSError as err:
        print '%s::Could not navigate to working directory.\n\t%s'\
                % (get_time(), err)
        sys.exit(-1)
    pfile = _parse_pfile(args.pf)
    print "%s::Checking parameter file." % get_time()
    _check_pfile(pfile)
    station_list = _create_station_list(args, pfile)
    input_files = ('mode_set.in',
                   'vgrids.in',
                   'interfaces.in',
                   'gridsave.in',
                   'frechet.in')
    for input_file in input_files:
        try:
            print '%s::Copying %s to working directory' \
                    % (get_time(), pfile[input_file])
            shutil.copyfile(pfile[input_file], input_file)
        except IOError as err:
            print "%s::Could not copy %s to working directory" \
                    % (get_time(), input_file)
            sys.exit(-1)
    print "%s::Creating dummy receivers.in file." % get_time()
    _create_dummy_receivers_file(pfile)
    _write_propgrid(pfile)
    gen_sta_tt_maps(station_list)
    for input_file in input_files:
        try:
            print '%s::Removing %s from working directory' \
                    % (get_time(), pfile[input_file])
            os.remove(input_file)
        except IOError as err:
            print "%s::Could not remove %s from working directory" \
                    % (get_time(), input_file)
    sys.exit(0)

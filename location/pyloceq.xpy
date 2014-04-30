def _main():
    """
    Standard main() function. Execution control begins here.
    """
    from anfseistools.ant import pf_2_cfg,\
                                     create_event_list,\
                                     write_origin
    from anfseistools.core import Locator,\
                                 parse_cfg
    from antelope.datascope import closing, dbopen
    args = _parse_command_line()
    pf_2_cfg(args.pf, 'pyloceq')
    cfg_dict = verify_config_file(parse_cfg('pyloceq.cfg'))
    locator = Locator(cfg_dict)
    with closing(dbopen(args.db, 'r+')) as db:
        tbl_event = db.schema_tables['event']
        if args.subset:
            view = tbl_event.join('origin')
            view = view.subset(args.subset)
            tbl_event = view.separate('event')
        for record in tbl_event.iter_record():
            evid = record.getv('evid')[0]
            view = tbl_event.subset('evid == %d' % evid)
            event_list = create_event_list(view)
            for event in event_list:
                origin = event.preferred_origin
                print('Relocating evid: %d'
                        % event.evid)
                origin = locator.locate_eq(origin)
                if origin == None:
                    print 'Could not relocate orid: %d' \
                            % event.preferred_origin.orid
                    continue
                origin.update_predarr_times(cfg_dict)
                write_origin(origin, db)
    return 0

def _parse_command_line():
    """
    Parse command line arguments. Return dictionary-like object
    containing results.
    """
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('db', type=str, help='Input/output databse.')
    parser.add_argument('-s', '--subset', type=str, help='Subset expression.')
    parser.add_argument('-p', '--pf', type=str, help='Parameter file.')
    return parser.parse_args()

if __name__ == '__main__': sys.exit(_main())
else: raise ImportError

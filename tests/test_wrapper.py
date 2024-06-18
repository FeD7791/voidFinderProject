from voidfindertk.zobov._wrapper import run_vozinit


def test_run_vozinit():

    input_params = {
        'position_file': 'tracers_zobov.raw',
        'buffer_size': '0.8',
        'box_size': '500',
        'number_of_divisions': '2',
        'suffix_describing_this_run': 'output'
    }
    len_box = '700070'
    generic_output_run_vozinit = [
        f'boz = {input_params["box_size"]}',
        f'np = {len_box}',
        f'np: {len_box}, x: 1.06e-06,0.999999; y: 8e-07,1; z: 8.4e-07,1',
        's = 0.0119048, bf = 0.8, g = 0.799911.',
        'b=(0,0,0), c=(0.25,0.25,0.25), nvp=85789, nvpbuf=700070',
        'b=(0,0,1), c=(0.25,0.25,0.75), nvp=82254, nvpbuf=700070',
        'b=(0,1,0), c=(0.25,0.75,0.25), nvp=92853, nvpbuf=700070',
        'b=(0,1,1), c=(0.25,0.75,0.75), nvp=83614, nvpbuf=700070',
        'b=(1,0,0), c=(0.75,0.25,0.25), nvp=88198, nvpbuf=700070',
        'b=(1,0,1), c=(0.75,0.25,0.75), nvp=84321, nvpbuf=700070',
        'b=(1,1,0), c=(0.75,0.75,0.25), nvp=96789, nvpbuf=700070',
        'b=(1,1,1), c=(0.75,0.75,0.75), nvp=86253, nvpbuf=700070',
        'Nvp range: 82254,96789',
        'Nvpbuf range: 700070,700070',
        'Writing script file to scroutput.',
        ''
    ]

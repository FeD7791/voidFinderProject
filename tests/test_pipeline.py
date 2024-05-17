from voidfindertk.data_box import DataBox
from voidfindertk.sphericalvf.spherical_tools import SphericalVF, spherical_void_finder, SphericalVoids


def test_sphericalvf_class(mkbox, make_spherical_voids_params):
    box = mkbox(seed=42, size=1000) #Build artificial box of size 1000
    params = make_spherical_voids_params() #Build params for spherical void output
    data_box = DataBox(box) #Transform box to universal Box
    spherical_void = SphericalVoids(**params) #Build artificial void output using params
    model_find_output = {'voids':spherical_void} #Artificial output for mk_vbox step
    model = SphericalVF()
    output = model.mk_vbox(model_find_output,data_box) #Build artificial mk_vbox
    assert box == data_box.box #For spherical, data_box is the same as box 
    assert data_box == model.preprocess(databox=data_box) #preprocess in spherical won't modify data_box 
    assert output['box'] == data_box.box 
    assert output['voids'] == spherical_void
    assert (len(box),len(spherical_void)) == output['sparse'].shape

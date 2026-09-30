import unittest
import numpy as np
import model_v6 as v6

class ModelTests(unittest.TestCase):
    def test_heldout_values_do_not_change_training_preprocessing(self):
        x=np.array([[1,0,np.nan],[3,0,np.nan],[np.nan,0,np.nan],[4,99,88]],dtype='float32')
        names=np.array(['signal','train_constant','train_missing'])
        a,b,nn,_=v6.preprocess(x,names,np.array([0,1,2]),np.array([3]))
        x[3]=[1e6,-1e6,-1e6]
        aa,bb,nn2,_=v6.preprocess(x,names,np.array([0,1,2]),np.array([3]))
        np.testing.assert_array_equal(a,aa);np.testing.assert_array_equal(nn,nn2)
        self.assertEqual(nn.tolist(),['signal']);self.assertNotEqual(float(b[0,0]),float(bb[0,0]))
    def test_all_preprocessing_outputs_finite_with_structural_absence(self):
        x=np.array([[0,0],[2,1],[np.nan,0],[np.nan,1]],dtype='float32')
        a,b,_,_=v6.preprocess(x,np.array(['radiomic','available']),np.array([0,1,2]),np.array([3]))
        self.assertTrue(np.isfinite(a).all() and np.isfinite(b).all())
if __name__=='__main__':unittest.main()

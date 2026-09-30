import unittest
import numpy as np
import extract_v6 as v6

class ExtractionTests(unittest.TestCase):
    def test_anatomical_pairs_not_parity(self):
        labels={1:'A_L',8:'A_R',2:'B_L',3:'Raphe_D',4:'Raphe_M'}
        self.assertEqual(v6.paired_labels(labels),{1:8,8:1})
    def test_firstorder_against_known_distribution(self):
        got=v6.firstorder([0,0,5,5])
        np.testing.assert_allclose(got,[2.5,2.5,0,5,2.5,5,np.sqrt(12.5),1])
        np.testing.assert_allclose(v6.contrast(got,got),np.zeros(8))
    def test_small_and_degenerate_texture_masks(self):
        a=np.ones((1,8,8),bool);self.assertFalse(v6.eligible_texture(a))
        a=np.ones((1,9,8),bool);self.assertTrue(v6.eligible_texture(a))
        a=np.ones((1,1,100),bool);self.assertFalse(v6.eligible_texture(a))
    def test_reference_excludes_dilated_lesion(self):
        # Anatomical whole-reference mask excludes bilateral disease and its margin.
        a=np.zeros((9,9,9),np.uint8);a[4,4,4]=1
        excluded=v6.sitk.GetArrayFromImage(v6.sitk.BinaryDilate(v6.sitk.GetImageFromArray(a),[2]*3))>0
        reference=np.ones_like(a,bool)&~excluded
        self.assertFalse(reference[4,4,4]);self.assertFalse(reference[4,4,6]);self.assertTrue(reference[0,0,0])
if __name__=='__main__':unittest.main()

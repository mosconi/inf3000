
import argparse
import numpy as np

import roadef
import instance
from instance import Instance,ValidateStatus

import unittest

class Instance_test(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser( description="Test instance",
                                   parents=[instance.parser,
                                            roadef.parser]
                                   )

        args = parser.parse_args(["-p","data/model_example.txt","-i","data/assignment_example.txt"])
        self.args = args

        self.inst = Instance(args)
        
    def tearDown(self):
        self.args.instance_filename.close()
        self.args.original_solution_filename.close()
        
    def test_map_validate_identity(self):
        x0 = self.inst.map_assign()
        v = self.inst.map_validate(x0)
        self.assertEqual(v.status,ValidateStatus.Valid)
        self.assertEqual(v.obj,4200)

    def test_map_validate_exemple1(self):
        x = np.zeros((self.inst.nproc,self.inst.nmach), dtype=np.int32)
        x[0,0] = 1
        x[1,2] = 1
        x[2,0] = 1

        v = self.inst.map_validate(x)
        self.assertEqual(v.status,ValidateStatus.Valid)
        self.assertEqual(v.obj,3510)

    def test_map_validate_exemple2(self):
        x = np.zeros((self.inst.nproc,self.inst.nmach), dtype=np.int32)
        x[0,0] = 1
        x[1,2] = 1
        x[2,1] = 1

        v = self.inst.map_validate(x)
        self.assertEqual(v.status,ValidateStatus.Valid)
        self.assertEqual(v.obj,2411)

    def test_map_validate_not_alloc(self):
        x = np.zeros((self.inst.nproc,self.inst.nmach), dtype=np.int32)

        v = self.inst.map_validate(x)
        self.assertEqual(v.status,ValidateStatus.Allocation)

    def test_map_validate_2_allocs(self):
        x = np.zeros((self.inst.nproc,self.inst.nmach), dtype=np.int32)
        x[0,0] =1
        x[0,1] =1
        x[1,3] =1
        x[2,0] =1
        v = self.inst.map_validate(x)
        self.assertEqual(v.status,ValidateStatus.Allocation)

    def test_map_validate_conflict(self):
        x = np.zeros((self.inst.nproc,self.inst.nmach), dtype=np.int32)
        x[0,0] =1
        x[1,0] =1
        x[2,0] =1
        v = self.inst.map_validate(x)
        self.assertEqual(v.status,ValidateStatus.Conflict)

    def test_map_validate_spread(self):
        x = np.zeros((self.inst.nproc,self.inst.nmach), dtype=np.int32)
        x[0,0] =1
        x[1,1] =1
        x[2,0] =1
        v = self.inst.map_validate(x)
        self.assertEqual(v.status,ValidateStatus.Spread)


if __name__ == '__main__':
    unittest.main()

import unittest 
from instance import Instance

model_example_str="""2
1 100
0 10
4
0 0 30 400 16 80 0 1 4 5
0 0 10 240 8 160 1 0 3 4
1 1 15 100 12 80 4 3 0 2
1 2 10 100 8 80 5 4 2 0
2
2 0
1 1 0
3
0 12 10 1000
0 10 20 100
1 6 200 1
1
0 1 20
10
1 10 100
"""
assign_example_str="0 3 0"
assign_example=[0, 3, 0]


class Test_Instance(unittest.TestCase):
    
    def setUp(self):
        self.inst = Instance(inline=True,
                             assign=assign_example_str,
                             model=model_example_str)

    def test_instance(self):
        self.assertIsInstance(self.inst,Instance)

    def test_assign(self):
        self.assertEqual(self.inst.assign(),[0,3,0])

    def test_validate0(self):
        self.assertTrue(self.inst.validate([0,3,0]))
        
    def test_objective0(self):
        self.assertEqual(self.inst.objective([0, 3, 0]),4200)

    def test_validate1(self):
        self.assertTrue(self.inst.validate([0,2,0]))

    def test_objective1(self):
        self.assertEqual(self.inst.objective([0, 2, 0]),3510)

    def test_validate2(self):
        self.assertTrue(self.inst.validate([0,2,1]))

    def test_objective2(self):
        self.assertEqual(self.inst.objective([0, 2, 1]),2411)

    def test_validate3(self):
        self.assertRaises(Exception,self.inst.validate,[0,0,0])

    def test_objective3(self):
        self.assertRaises(Exception,self.inst.objective,[0, 0, 0])

    def test_validate4(self):
        self.assertRaises(Exception,self.inst.validate,[2,3,0])

if __name__ == '__main__':
    unittest.main()

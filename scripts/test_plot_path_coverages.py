import unittest, importlib
from tempfile import TemporaryDirectory
from collections import defaultdict
import io
from contextlib import redirect_stdout
from plot_path_coverages import process_bubble, add_stats


class Test_process_bubble(unittest.TestCase):
	def test_process_bubble1(self):
		lines = ["path0\t0.1\t0\t31", "path1\t0\t1\t31" ]
		path_to_fraction, consider_variant = process_bubble(lines)
		self.assertTrue("path0" in path_to_fraction)
		self.assertTrue("path1" in path_to_fraction)
		self.assertTrue(path_to_fraction["path0"] == 0.1)
		self.assertTrue(path_to_fraction["path1"] == 0)
		self.assertFalse(consider_variant)
	def test_process_bubble2(self):
		lines = ["path0\t-1\t0\t31", "path1\t-1\t1\t31" ]
		path_to_fraction, consider_variant = process_bubble(lines)
		self.assertTrue("path0" in path_to_fraction)
		self.assertTrue("path1" in path_to_fraction)
		self.assertTrue(path_to_fraction["path0"] == -1)
		self.assertTrue(path_to_fraction["path1"] == -1)
		self.assertFalse(consider_variant)
	def test_process_bubble3(self):
		lines = [
				"path0\t0.95\t0\t31",
				"path1\t0.95\t0\t31",
				"path2\t0.95\t0\t31",
				"path3\t0.95\t0\t31",
				"path4\t0.95\t0\t31",
				"path5\t0.95\t0\t31",
				"path6\t0.95\t0\t31",
				"path7\t0.95\t0\t31",
				"path8\t0.95\t0\t31",
				"path9\t0.95\t0\t31",
				"path10\t0.95\t0\t31",
				"path11\t0.95\t0\t31",
				"path12\t0.95\t0\t31",
				"path13\t0.95\t0\t31",
				"path14\t0.95\t0\t31",
				"path15\t0.95\t0\t31",
				"path16\t0.95\t0\t31",
				"path17\t0.95\t0\t31",
				"path18\t0.95\t0\t31",
				"path19\t0.60\t1\t31",
				 ]
		path_to_fraction, consider_variant = process_bubble(lines)
		for i in range(0,19):
			self.assertTrue("path" + str(i) in path_to_fraction)
			self.assertTrue(path_to_fraction["path" + str(i)] == 0.95)
		self.assertTrue("path" + str(i) in path_to_fraction)
		self.assertTrue(path_to_fraction["path19"] == 0.60)
		self.assertTrue(consider_variant)

	def test_process_bubble4(self):
		lines = [
				"path0\t0.95\t0\t31",
				"path1\t0.95\t0\t31",
				"path2\t0.95\t0\t31",
				"path3\t0.95\t0\t31",
				"path4\t0.95\t0\t31",
				"path5\t0.95\t0\t31",
				"path6\t0.95\t0\t31",
				"path7\t0.95\t0\t31",
				"path8\t0.95\t0\t31",
				"path9\t0.95\t0\t31",
				"path10\t0.95\t0\t31",
				"path11\t0.95\t0\t31",
				"path12\t0.95\t0\t31",
				"path13\t0.95\t0\t31",
				"path14\t0.95\t0\t31",
				"path15\t0.95\t0\t31",
				"path16\t0.95\t0\t31",
				"path17\t0.95\t0\t31",
				"path18\t0.95\t0\t31",
				"path19\t0.40\t1\t31",
				 ]
		path_to_fraction, consider_variant = process_bubble(lines)
		for i in range(0,19):
			self.assertTrue("path" + str(i) in path_to_fraction)
			self.assertTrue(path_to_fraction["path" + str(i)] == 0.95)
		self.assertTrue("path" + str(i) in path_to_fraction)
		self.assertTrue(path_to_fraction["path19"] == 0.40)
		# allele is rare but not well covered
		self.assertFalse(consider_variant)


class Test_add_stats(unittest.TestCase):
	def test_add_stats1(self):
		path_to_fraction_all = defaultdict(list)
		path_to_fraction_present = defaultdict(list)
		lines = [
				"path0\t0.95\t0\t31",
				"path1\t0.95\t0\t31",
				"path2\t0.95\t0\t31",
				"path3\t0.95\t0\t31",
				"path4\t0.95\t0\t31",
				"path5\t0.95\t0\t31",
				"path6\t0.95\t0\t31",
				"path7\t0.95\t0\t31",
				"path8\t0.95\t0\t31",
				"path9\t0.95\t0\t31",
				"path10\t0.95\t0\t31",
				"path11\t0.95\t0\t31",
				"path12\t0.95\t0\t31",
				"path13\t0.95\t0\t31",
				"path14\t0.95\t0\t31",
				"path15\t0.95\t0\t31",
				"path16\t0.95\t0\t31",
				"path17\t0.95\t0\t31",
				"path18\t0.95\t0\t31",
				"path19\t0.60\t1\t31",
				 ]
		add_stats(lines, path_to_fraction_all, path_to_fraction_present)

		for i in range(0,19):
			self.assertTrue("path" + str(i) in path_to_fraction_all)
			self.assertTrue(path_to_fraction_all["path" + str(i)] == [0.95])
		self.assertTrue("path" + str(i) in path_to_fraction_all)
		self.assertTrue(path_to_fraction_all["path19"] == [0.60])

		for i in range(0,19):
			self.assertTrue("path" + str(i) in path_to_fraction_present)
			self.assertTrue(path_to_fraction_present["path" + str(i)] == [0.95])
		self.assertTrue("path" + str(i) in path_to_fraction_present)
		self.assertTrue(path_to_fraction_present["path19"] == [0.60])


	def test_add_stats2(self):
		path_to_fraction_all = defaultdict(list)
		path_to_fraction_present = defaultdict(list)
		lines = ["path0\t-1\t0\t31", "path1\t-1\t1\t31" ]
		add_stats(lines, path_to_fraction_all, path_to_fraction_present)

		self.assertTrue("path0" in path_to_fraction_all)
		self.assertTrue("path1" in path_to_fraction_all)
		self.assertTrue(path_to_fraction_all["path0"] == [-1])
		self.assertTrue(path_to_fraction_all["path1"] == [-1])

		self.assertFalse(path_to_fraction_present)
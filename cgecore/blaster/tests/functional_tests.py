import unittest
from subprocess import PIPE, run
import sys
import os.path
import shutil

# This is not best practice but for testing, this is the best I could come up
# with
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from blaster import Blaster


python3 = "python3"
root_test_dir = os.path.abspath(os.path.dirname(__file__))
run_test_dir = root_test_dir + "/running_test"

test_names = ["test1", "test1_db", "test2", "test2_db"]
test_data = {
    # Test hit to test_db01
    test_names[0]: root_test_dir + "/data/test_isolate_01.fa",
    test_names[1]: root_test_dir + "/data/test_db_01.fa",
    test_names[2]: root_test_dir + "/data/test_isolate_01.fa",
    test_names[3]: root_test_dir + "/data/test_db_01.fa",
}
working_dir = os.path.dirname(os.path.realpath(__file__))


class ResFinderRunTest(unittest.TestCase):

    def setUp(self):
        # Change working dir to test dir
        os.chdir(working_dir)
        # Does not allow running two tests in parallel
        os.makedirs(run_test_dir, exist_ok=False)

    def tearDown(self):
        shutil.rmtree(run_test_dir)

    def test_run_blaster_script(self):
        test_dir = run_test_dir + "/" + test_names[0]
        os.makedirs(test_dir)

        inputfile = test_data[test_names[0]]

        dirs, filename = os.path.split(test_data[test_names[1]])
        db_name, ext = os.path.splitext(filename)
        databases = [db_name]
        db_path = dirs
        out_path = test_dir
        min_cov = 0.6
        threshold = 0.8
        blast = "blastn"

        cmd = (python3 + " ../run_blaster.py"
               + " --input " + inputfile
               + " --databases " + db_name
               + " --out_path " + out_path
               + " --blastPath " + blast
               + " --databasePath " + db_path
               + " --min_cov 0.6"
               + " --threshold 0.8")

        procs = run(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        # Expected output files

        in_dirs, in_filename = os.path.split(test_data[test_names[0]])
        in_name, ext = os.path.splitext(in_filename)

        out_aln = test_dir + "/" + in_name + "_hit_alignments.txt"
        out_res = test_dir + "/" + in_name + "_results.txt"
        out_xml = test_dir + "/tmp/out_" + db_name + ".xml"

        with open(out_aln, "r") as fh:
            fh.readline()
            fh.readline()
            check_result = fh.readline()
        self.assertIn("blaB-2_1_AF189300", check_result)
        with open(out_res, "r") as fh:
            fh.readline()
            fh.readline()
            fh.readline()
            fh.readline()
            check_result = fh.readline()
        self.assertIn("blaB-2_1_AF189300", check_result)
        with open(out_xml, "r") as fh:
            for i in range(0, 9):
                check_result = fh.readline()
        self.assertIn("blaB-2_1_AF189300", check_result)

    def test_do_a_simple_blast(self):
        # Maria has an E. coli isolate, with beta lactam resistance.
        # She wants to blast the isolate against her own beta lactam database.

        test_dir = run_test_dir + "/" + test_names[2]
        os.makedirs(test_dir)

        # Maria sets up some variables in her script.
        inputfile = test_data[test_names[2]]

        dirs, filename = os.path.split(test_data[test_names[3]])
        db_name, ext = os.path.splitext(filename)
        databases = [db_name]
        db_path = dirs
        out_path = test_dir
        min_cov = 0.6
        threshold = 0.8
        blast = "blastn"

        # Maria uses the blaster method from the module to do her BLASTing.
        blast_run = Blaster(inputfile=inputfile, databases=databases,
                            db_path=db_path, out_path=out_path,
                            min_cov=min_cov, threshold=threshold, blast=blast)

        results = blast_run.results
        query_align = blast_run.gene_align_query
        homo_align = blast_run.gene_align_homo
        sbjct_align = blast_run.gene_align_sbjct

        tab_file = "%s/%s_results.txt" % (out_path, filename)
        tab = open(tab_file, "w")
        for db in results:
          tab.write("%s\n" % (db))
          tab.write("Hit\tIdentity\tAlignment Length/Gene Length\tPosition in"
                    " reference\tContig\tPosition in contig\n")

          for hit in results[db]:
             header = results[db][hit]["sbjct_header"]
             ID = results[db][hit]["perc_ident"]
             sbjt_length = results[db][hit]["sbjct_length"]
             HSP = results[db][hit]["HSP_length"]
             positions_contig = "%s..%s" % (results[db][hit]["query_start"],
                                            results[db][hit]["query_end"])
             positions_ref = "%s..%s" % (results[db][hit]["sbjct_start"],
                                         results[db][hit]["sbjct_end"])
             contig_name = results[db][hit]["contig_name"]

             self.assertIn("blaB-2_1_AF189300", header)

             # Write tabels
             tab.write("%s\t%.2f\t%s/%s\t%s\t%s\t%s\n" % (header, ID, HSP,
                                                          sbjt_length,
                                                          positions_ref,
                                                          contig_name,
                                                          positions_contig))

        tab.close()


if __name__ == "__main__":
    unittest.main()

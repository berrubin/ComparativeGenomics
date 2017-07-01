from Bio.Phylo.PAML import yn00
from Bio import SeqIO
from Bio.Phylo.PAML import codeml
import subprocess
import os

def branch_site_d_worker(orthogroup, workingdir):
    print orthogroup
    template = open("/Genomics/kocherlab/berubin/annotation/orthology/codeml.CladeD.ctl", 'rU')
    os.chdir(workingdir)
    if not os.path.isdir("og_%s_working" % orthogroup):
        os.mkdir("og_%s_working" % orthogroup)
    os.chdir("og_%s_working" % orthogroup)
    alt_ctl = open("og_%s_alt.ctl" % (orthogroup), 'w')
    for line in template:
        if "seqfile = " in line:
            alt_ctl.write(line.replace("dummy_file", "../og_cds_%s.afa" % (orthogroup)))
        elif "treefile = " in line:
            alt_ctl.write(line.replace("dummy_file", "../og_%s.tree" % (orthogroup)))
        elif "outfile = " in line:
            alt_ctl.write(line.replace("dummy_file", "../og_%s.alt" % (orthogroup)))
        else:
            alt_ctl.write(line)
    template.close()
    alt_ctl.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", "og_%s_alt.ctl" % (orthogroup)]
    subprocess.call(cmd)
    alt_ctl = open("og_%s_alt.ctl" % ( orthogroup), 'rU')
    nul_ctl = open("og_%s_nul.ctl" % ( orthogroup), 'w')
    for line in alt_ctl:
        if "model = " in line:
            nul_ctl.write(line.replace("3", "0"))
        elif "outfile = " in line:
            nul_ctl.write(line.replace("alt", "nul"))
        else:
            nul_ctl.write(line)
    alt_ctl.close()
    nul_ctl.close()
    cmd = ["/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", "og_%s_nul.ctl" % ( orthogroup)]
    subprocess.call(cmd)
    os.chdir("/Genomics/kocherlab/berubin/annotation/orthology")
    
#this tests for selection to be positive on foreground branches
def branch_positive_worker(orthogroup, workingdir):
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.alt" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=2, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 0, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = False)
    cml.set_options(NSsites=[0])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.nul" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=2, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 1, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6)
    cml.set_options(NSsites=[0])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)

#this just tests that foreground branches are different than overall branches
def branch_worker(orthogroup, workingdir):
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.alt" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=2, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 0, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = False)
    cml.set_options(NSsites=[0])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.nul" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=0, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 0, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = False)
    cml.set_options(NSsites=[0])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)

def branch_site_worker(orthogroup, workingdir):
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.alt" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=2, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 0, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = False)
    cml.set_options(NSsites=[2])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.nul" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=2, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 1, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6)
    cml.set_options(NSsites=[2])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)

def free_ratios_worker(orthogroup, workingdir):
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.alt" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=1, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 0, omega = 1, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = False)
    cml.set_options(NSsites=[0])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)


def ancestor_reconstruction(orthogroup, workingdir):
    cml = codeml.Codeml(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), tree = "%s/og_%s.tree" % (workingdir, orthogroup), out_file = "%s/og_%s.anc" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    cml.set_options(runmode=0,fix_blength=0,seqtype=1,CodonFreq=2, model=0, icode=0, clock = 0, aaDist=0, Mgene = 0, fix_kappa = 0, kappa = 2, fix_omega = 0, omega = 1, getSE = 0, RateAncestor = 1, cleandata = 0, Small_Diff = .45e-6, verbose = False)
    cml.set_options(NSsites=[0])
    cml.print_options()
    cml.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/codeml", verbose = False)


def pairwise_yn(orthogroup, workingdir):
    print workingdir
    yn = yn00.Yn00(alignment = "%s/og_cds_%s.afa" % (workingdir, orthogroup), out_file = "%s/og_%s.yn" % (workingdir, orthogroup), working_dir = "%s/og_%s_working" % (workingdir, orthogroup))
    yn.set_options(icode = 0, verbose = 0, weighting = 0, commonf3x4 = 0)
    yn_results = yn.run(command = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin/yn00", verbose = False)
    return yn_results

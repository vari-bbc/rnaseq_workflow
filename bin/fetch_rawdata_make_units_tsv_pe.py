import argparse
import glob
import os
import re
import shutil
from typing import List
from typing import Dict


def write_units_tsv(sample_name_dict: Dict) -> None:
    fout_name = os.path.join("bin", "units.tsv")
    shutil.copy(fout_name, os.path.join("bin", "units.tsv.bak"))

    with open(fout_name, "w") as fout:
        fout.write(f"sample\tfq1\tfq2\n")
        for sample in sample_name_dict:
            fout.write(
                f"{sample}\t{sample_name_dict[sample][0]}\t{sample_name_dict[sample][1]}\n"
            )
    print("writing units.tsv is done!")


def get_sample_name_list(fq_list: List) -> Dict:
    # key sample_name : value [fq1, fq2]
    results_dict = dict()
    for fq in fq_list:
        fq1: str = ""
        fq2: str = ""
        sample_name: str = re.compile("^(\w+)_L00\d").match(fq).group(1)
        if "_R1_" in fq and not "Undetermined" in fq:
            fq1: str = fq
            if not sample_name in results_dict:
                results_dict[sample_name] = [fq1, False]
            else:
                if results_dict[sample_name][0] is False:
                    results_dict[sample_name][0] = fq1
        elif "_R2_" in fq and not "Undetermined" in fq:
            fq2: str = fq
            if not sample_name in results_dict:
                results_dict[sample_name] = [False, fq2]
            else:
                if results_dict[sample_name][1] is False:
                    results_dict[sample_name][1] = fq2
        else:
            raise ValueError(f"{fq} doesn't have _R1_ or _R2_")

    return results_dict


def get_fq_list_from_genomics_core(fq_path: str) -> List[str]:
    """
    feed abs path for this function
    :param fq_path:
    :return: fq files list
    """
    cur_dir: str = os.getcwd()
    print(f"moving to {fq_path}")

    os.chdir(fq_path)
    fq_list: List[str] = sorted(glob.glob("*fastq.gz"))
    for fq in fq_list:
        if "Undetermined" in fq:
            fq_list.remove(fq)

    print(f"moving back to {cur_dir}")
    os.chdir(cur_dir)
    return fq_list


def copy_fq_to_raw_dir(fq_path: str, fq_list: List, dest_dir: str) -> None:
    # dest_dir: str = "raw_data/"
    assert os.path.isdir(dest_dir)

    for fq in fq_list:
        full_path: str = os.path.join(fq_path, fq)
        assert os.path.isfile(full_path)

        print(f"copying {full_path} -> {dest_dir}")
        shutil.copy(full_path, dest_dir)
    print(f"copying fq files is done.")


def main(fq_path: str, dest_dir: str) -> None:
    fq_list = get_fq_list_from_genomics_core(fq_path=fq_path)
    copy_fq_to_raw_dir(fq_path=fq_path, fq_list=fq_list, dest_dir=dest_dir)
    sample_name_list = get_sample_name_list(fq_list=fq_list)
    write_units_tsv(sample_name_dict=sample_name_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Python script for fetching fastq.gz files and create units.tsv"
    )

    parser.add_argument(
        "-i",
        "--fq_dir",
        action="store",
        type=str,
        required=True,
        help="fastq.gz dir",
    )

    parser.add_argument(
        "-d",
        "--dest_dir",
        action="store",
        type=str,
        default="raw_data",
        help="destination dir",
    )

    args = parser.parse_args()
    main(args.fq_dir, args.dest_dir)

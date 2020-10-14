#!/usr/bin/python

import sys, getopt
import threading
import pathlib
from time import sleep
from time import time
from argparse import ArgumentParser
import io
import subprocess

import pandas as pd
import logging
import pycurl
import statistics


output_directory = 'data/'
#base_url = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/'
#https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR8296149&format=fasta
#http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRRNNNNNN
#https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/SRR8296149/SRR8296149.1


finished_file_bytes = 0
thread_list = []
concurrency = 1


lock_sample_list = threading.Lock()
lock_active_transfer_list = threading.Lock()


kb = 1024
# callback function for c.XFERINFOFUNCTION
def status(download_t, download_d, upload_t, upload_d):
    #STREAM.write
    print('Downloading: {}/{} kiB ({}%)\r'.format(
        str(int(download_d/kb)),
        str(int(download_t/kb)),
        str(int(download_d/download_t*100) if download_t > 0 else 0)
    ))
    #STREAM.flush()


def download_monitor (target_throughput, sample_list, active_transfer_list):
    global finished_file_bytes
    last_size = 0
    last_concurrency_update_t = 0
    throughput_list = []
    while True:
        sleep(1)
        size = finished_file_bytes
        if len(sample_list) == 0 and len(active_transfer_list) == 0:
           return
        for f in active_transfer_list:
            try:
                size += pathlib.Path(f).stat().st_size
            except:
                continue
        throughput = ((size-last_size)*8)/(1000*1000.0)
        if throughput == 0:
            continue
        throughput_list.append(throughput)
        print ("Throughput {} Mbps, target {} Mbps, Active files: {}, Remaining files {}".format(throughput, target_throughput, len(active_transfer_list), len(sample_list)))
        last_size = size
        if len(throughput_list) > 5 and statistics.mean(throughput_list[-5:]) < 0.9* target_throughput \
                and time() > last_concurrency_update_t + 6:
            if len(sample_list) < 2:
                continue
            print("Creating new thread " + str(statistics.mean(throughput_list[-5:])))
            thread = threading.Thread(target=file_downloader, args=(sample_list, active_transfer_list), daemon=True)
            thread.start()
            thread_list.append(thread)
            last_concurrency_update_t = time()
    return



def file_downloader (sample_list, active_transfer_list):
    global finished_file_bytes
    print("Running thread...")
    while True:
        lock_sample_list.acquire()
        if len(sample_list) == 0:
            print("Exiting thread...")
            return

        curl = pycurl.Curl()
        # Fetch a file from file list and add is to currently transferred file list
        # Synchronized operation due to using concurrency

        filename = sample_list.pop()
        lock_sample_list.release()

        lock_active_transfer_list.acquire()
        file_path = output_directory + filename
        active_transfer_list.append(file_path)
        lock_active_transfer_list.release()




        # initialize and start curl file download
        sample_ftp_url = discover_ftp_paths([filename])[0]
        #print("Starting to download " + filename, sample_ftp_url)
        fp = open(file_path, "wb")
        curl.setopt(pycurl.URL, sample_ftp_url)
        curl.setopt(pycurl.WRITEDATA, fp)
        #curl.setopt(curl.NOPROGRESS, False)
        #curl.setopt(curl.XFERINFOFUNCTION, status)
        try:
            curl.perform()
        except pycurl.error as exc:
            raise ValueError("Unable to reach %s (%s)" % (sample_ftp_url, exc))
        file_size = pathlib.Path(file_path).stat().st_size
        print("Finished {} size: {} MB".format(filename, (file_size/(1024*1024))))
        curl.close()
        fp.close()

        # Remove file from currently transferred file list
        lock_active_transfer_list.acquire()
        active_transfer_list.remove(file_path)
        finished_file_bytes += file_size
        lock_active_transfer_list.release()

#https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&fields=fastq_ftp&accession=
def discover_ftp_paths (sample_list):

    command = "srapath " + " ".join(sample_list)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    paths = result.stdout.decode('utf-8').splitlines()
    #print ("command", command, " paths:", paths, "hah")
    return paths


fileList = []

if __name__ == "__main__":

    parser = ArgumentParser()

    parser.add_argument('-i', '--input',
                  action="store", dest="sample_list",
                  help="input file for sample list", default="samples.tsv")
    parser.add_argument('-o', '--output',
                      action="store", dest="output_directory",
                      help="output directory to save sample files", default="data")
    parser.add_argument('-c', '--concurrency',
                      action="store", dest="concurrency", type=int,
                      help="default concurrency value", default=1)
    parser.add_argument('-t', '--target',
                      action="store", dest="target_throughput", type=int,
                      help="target throughput for the transfer", default=0)

    args = parser.parse_args()

    sample_file = pd.read_csv(args.sample_list, sep="\s+", dtype=str).set_index("sample", drop=False)
    sample_list = sample_file["sample"].values.tolist()

    #sample_list_ftp_url = discover_ftp_paths (sample_list)



    #sys.exit(-1)
    active_transfer_list = []

    pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)
    t = threading.Thread(target=download_monitor, \
                     args=(args.target_throughput, sample_list, active_transfer_list,), daemon=True)
    t.start()
    for i in range (args.concurrency):
        thread = threading.Thread(target=file_downloader, \
                                  args=(sample_list, active_transfer_list,), daemon=True)
        logging.info("Main    : starting thread " + str(i))
        thread.start()
        thread_list.append(thread)
    t.join()
import re
import subprocess
from dataclasses import dataclass

from adapters.primer3_to_ipcress import Primer3ToIpcressAdapter
from utils.write_output_files import write_ipcress_input

@dataclass
class IpcressParams:
    pass

@dataclass
class IpcressResult:
    stnd: bytes
    err: bytes


class Ipcress:
    def __init__(self, params):
        self.params = params

    def run(self):
        print("Running iPCRess...")
        if self.params["verbose"]:
            print(self.params["verbose"])
            print('iPCRess params:')
            print(self.params)
        params = self.params

        if params['primers']:                               # Comes from --primers argument, should be a path to a txt file
            if params['verbose']:
                print('Loading custom iPCRess input file')
            input_path = params['primers']
            result = self.run_ipcress(input_path, params)
        else:
            if params['verbose']:
                print('Building iPCRess input file.')

            adapter = Primer3ToIpcressAdapter()
            adapter.prepare_input(
                params['p3_csv'], params['min'], params['max'], params['dir']
            )
            input_path = write_ipcress_input(params['dir'], adapter.formatted_primers)

            result = self.run_ipcress(input_path, params)
            if params['verbose']:
                print(result.stnd.decode())
            self.validate_primers(
                result.stnd.decode(), adapter.primer_data, params
            )

        if params['verbose']:
            print('Finished!')

        return result

    def run_ipcress(self, input_path, params) -> IpcressResult:
        cmd = ' '.join([
            'ipcress', input_path, params['fasta'],'--mismatch', params['mismatch']
        ])

        cmd = self.prettify_output(params['pretty'], cmd)

        if params['verbose']:
            print('Running Exonerate iPCRess with the following command:\n{0}'.format(cmd))

        # ipcress = subprocess.Popen(
        #     cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        # )

        # stnd, err = ipcress.communicate()
        result = subprocess.run(cmd, shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        stnd = result.stdout
        err = result.stderr
        if not stnd:
            raise Exception(err)
        
        return IpcressResult(stnd, err)

    @staticmethod
    def validate_primers(ipcress_output, primer_data, params):
        if params["pretty"]:
            print('Output is pretty, skipping validation')
            return
        if params['verbose']:
            print('Validating primers...')

        for primer_pair in primer_data.keys():
            fwd_coord = primer_data[primer_pair]['F']['start']
            rev_coord = primer_data[primer_pair]['R']['start']

            match = re.search((
                fr'{primer_pair} \d+ A {fwd_coord} 0 '
                fr'B {rev_coord} 0 forward$'), ipcress_output)

            if not match:
                print(f'No valid primer pair found for {primer_pair}')


    @staticmethod
    def prettify_output(prettify, cmd):
        pretty_opt = 'false'
        if prettify:
            pretty_opt = 'true'

        return cmd + ' --pretty ' + pretty_opt

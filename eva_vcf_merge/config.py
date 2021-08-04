
class MergeConfig:
    def __init__(self, nextflow_binary, nextflow_config, bcftools_binary, output_dir):
        self.nextflow_binary = nextflow_binary
        self.nextflow_config = nextflow_config
        self.bcftools_binary = bcftools_binary
        self.output_dir = output_dir

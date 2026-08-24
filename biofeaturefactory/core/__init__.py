# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""BioFeatureFactory core: the input-producing stages.

Every module here is a CLI that PRODUCES the artifacts the feature pipelines
consume -- coordinate mappings and sequence records (variant_mapping), VCFs
(vcf_converter), and protein/codon alignments (msa_generation_pipeline,
codon_msa_pipeline). They are orchestrated by CLI, including from Nextflow, not
imported by the pipelines.

Layering: lib -> core -> pipelines. Adding a feature pipeline means writing one
consumer against the artifacts core already emits; it requires no change here
and none in lib.
"""

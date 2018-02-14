.. _data_organization:

=================
Data organisation
=================

Once all tasks have completed, the data will be uploaded into a
database. Currently, both sqlite_ and postgresql_ have been tested,
but mysql_ should work in principle as well.

The data is organised into a simple collection of tables. A benchmark
is composed of a single :term:`run` that groups together several
instances. Each :term:`instance` is a combination of a :term:`metric`,
:term:`tool`, input data and options. Each metric outputs a
tab-separated table that is uploaded into a separate table and adds
runtime execution performance into tables called :term:`_timings`.  A
typical collection of tables looks like this::

     > .tables
     # Maintenance tables
     run
     instance                     
     # Timing tables
     metric_timings
     tool_timings                             
     # metric tables
     bedtools_stats_allele_frequency
     bedtools_jaccard                    
     bcftools_stats_depth_distribution
     bcftools_stats_indel_context_length
     bcftools_stats_indel_context_summary
     bcftools_stats_indel_distribution
     bcftools_stats_quality
     bcftools_stats_singleton_stats
     bcftools_stats_substitution_types
     bcftools_stats_summary_numbers
     vcftools_tstv_by_count              
     vcftools_tstv_summary               

Table overview
===================

.. glossary::

   run
       Information about a benchmark run. Columns:

       id 
          Identification number of this run
       author
          The user name of the person running the pipeline
       created
          Date the benchmark run was created
       pipeline_name
          The name of the pipeline
       pipeline_version
          The pipeline version (git commit), typically
          the current git commit.
       config
          The benchmark configuration file in json format
       title
          The title of the benchmark run, see :ref:`configuration`
       description
          The description of the benchmark run, see :ref:`configuration`

   instance
       An instantiation of a combination of a particular
       :term:`metric`, :term:`tool`, input data and options.
	
       id
           Identification number of this instance
       run_id
           Reference to run
       completed
           Time that computation was completed
       input
           Input data
       metric_name
           Name of the metric
       metric_version
           Version of the metric
       metric_options
           Options supplied to the metric
       tool_name
           Name of the tool
       tool_version
           Tool version
       tool_options
           Options supplied to the tool
       meta_data
           Other environment variables

   timings
      Timing information

      instance_id
           Reference to instance
       host
           Execution host
       started
           Time that job was submitted
       completed
           Time that job was completed
       total_t
           Total time of job, including waiting in the queue
       wall_t
           Time spend in user/system in total
       user_t
           Time spend in user in job script
       sys_t
           Time spend in system in job script
       child_user_t
           Time spend in user in child processes. This is typically
           the tool/metric being executed
       child_sys_t
           Time spend in system in child processes. This is typically
           the tool/metric being executed
       statement
           Command line statement

   tags
      List of tags

      run_id
         Reference to the run
      tag
         A tag associated with the run

   arvados_job
      Arvados job information. This table is only present if arvados
      is ``--engine=arvados`` has been used

      run_id
         Reference to the run
      owner_uuid
         Arvados :term:`UUID` of the owner
      job_uuid
         Arvados :term:`UUID` of the job
      output_uuid
         Arvados :term:`UUID` of the output

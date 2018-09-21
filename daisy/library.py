"""Context manager for Library processing
=========================================

API
---
"""

import shutil
import os
import glob
import re

import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import cgatcore.pipeline as P
# import CGATCore.Arvados as Arvados

try:
    import arvados
    HAS_ARVADOS = True
except ImportError:
    HAS_ARVADOS = False


class LibraryContext(object):

    def __init__(self, params, options, args, argv, framework):
        self.params = params
        self.options = options
        self.args = args
        self.argv = argv
        self.workingdir_is_local = False
        self.shared_workingdir = None
        self.framework = framework
        self.logger = P.get_logger()

    def __exit__(self, tpe, value, tb):

        # TODO: handling exceptions in pipeline creation
        if tpe:
            return False

        if self.options.engine == "arvados":
            retval = True
        else:
            self.logger.debug("starting workflow")
            retval = P.main(self.options, self.args)

        if self.workingdir_is_local:
            # moving files to local directory - assumes that pipeline.log is
            # not in shared_workingdir
            self.logger.debug("moving files from {} to {}".format(
                self.shared_workingdir, self.workingdir))
            for root, dirs, files in os.walk(self.shared_workingdir):
                destdir = os.path.join(self.workingdir, os.path.relpath(
                    root, self.shared_workingdir))
                for d in dirs:
                    if not os.path.exists(os.path.join(destdir, d)):
                        os.makedirs(os.path.join(destdir, d))
                for fn in files:
                    shutil.move(os.path.join(root, fn),
                                os.path.join(destdir, fn))

            try:
                shutil.rmtree(self.shared_workingdir)
            except OSError as ex:
                self.logger.warn("could not remove {}: {}".format(self.shared_workingdir, ex))

        return retval

    def __enter__(self):

        options = self.options
        argv = self.argv
        params = self.params
        args = self.args
        framework = self.framework

        # special execution processes
        if options.engine == "arvados":
            # execute workflow as a crunch-script
            #
            # add explicit path to config file. In the future, upload
            # config file to arvados?
            if "--config" in argv or "--config-file" in argv:
                # replace config file name with absolute path
                for index, item in enumerate(argv):
                    if ".yml" in item:
                        try:
                            os.path.exists(item)
                        except:
                            pass
                        else:
                            argv[index] = os.path.abspath(item)
            else:
                assert "config-file" not in "".join(argv)
                argv.insert(1, "--config-file={}".format(
                    os.path.abspath(options.config_file)))

            crunch_json = Arvados.build_crunch_script(
                argv,
                yaml_file=options.config_file,
                project_uuid=options.project_uuid,
                framework=framework)

            statement = 'arv-crunch-job --job="$(cat {})"'.format(crunch_json)
            self.logger.debug("executing: {}".format(statement))
            retval = E.run(statement, return_stderr=True)

            # get crunch output
            try:
                m = re.search("job output (\w+\+\w+)", retval)
            except AttributeError:
                self.logger.warn("job output for Arvados crunch job is not found")
                output_hash = "unknown"
            else:
                output_hash = m.group(1)
            self.logger.debug("crunch job output hash: {}".format(output_hash))

            try:
                m = re.search("log collection is (\w+\+\w+)", retval)
            except AttributeError:
                self.logger.warn("log collection for Arvados crunch job is not found")
                log_hash = "unknown"
            else:
                log_hash = m.group(1)
            self.logger.debug("cruch job log hash: {}".format(log_hash))

            os.unlink(crunch_json)

            # move output and log to project specified (temporary solution)
            self.logger.info("move output to project {}".format(options.project_uuid))
            if output_hash != "unknown":
                api = arvados.api()
                output_uuid = Arvados.portable_data_hash2uuid(output_hash, api)
                api.collections().update(uuid=output_uuid,
                                         body={
                                             "collection": {
                                                 "name": "{}_output".format(options.job_name),
                                                 "owner_uuid": options.project_uuid
                                             }
                                         }).execute()
            self.logger.info("move log to project {}".format(options.project_uuid))
            if log_hash != "unknown":
                api = arvados.api()
                log_uuid = Arvados.portable_data_hash2uuid(log_hash, api)
                api.collections().update(uuid=log_uuid,
                                         body={
                                             "collection": {
                                                 "name": "{}_log".format(options.job_name),
                                                 "owner_uuid": options.project_uuid
                                             }
                                         }).execute()

        else:

            # sort out working directory
            if options.work_dir is not None:
                if not os.path.exists(options.work_dir):
                    shutil.makedirs(options.work_dir)
                self.original_dir = os.getcwd()
                os.chdir(options.work_dir)
                P.PARAMS["workingdir"] = options.work_dir

            self.workingdir = os.path.abspath(P.PARAMS["workingdir"])
            self.workingdir_is_local = (IOTools.is_local(self.workingdir) and
                                        P.will_run_on_cluster(P.PARAMS) and
                                        not options.without_cluster)

            # Special case: we want to upload with crunch
            if "keep" in args or "keep-and-load" in args:
                # check that this is a crunch script
                try:
                    job_info = arvados.current_job()
                except KeyError:
                    raise ValueError("keep called, but not a crunch context, "
                                     "use --engine=arvados")

                # link/mv files from dir to current dir
                # derive path from config_file set explicitely above.
                src_path = os.path.dirname(options.config_file)
                dest_path = os.path.abspath(P.PARAMS["workingdir"])
                files = ["pipeline.log", os.path.basename(options.config_file)]
                files.extend([os.path.basename(x)
                              for x in glob.glob(os.path.join(src_path, "*.dir"))])
                self.logger.debug("linking {} files/directories from {} to working directory {}".format(
                    len(files), src_path, dest_path))
                for fn in files:
                    src = os.path.join(src_path, fn)
                    if os.path.exists(src) and os.path.isfile(src):
                        os.symlink(src,
                                   os.path.join(dest_path, fn))
                    else:
                        # arvados crunch does not follow symlinks for directories, so
                        # copy.
                        shutil.copytree(src, os.path.join(dest_path, fn))

                # do not move to a shared local directory and turn cluster off
                self.workingdir_is_local = False
                P.PARAMS["without_cluster"] = True
                if "keep" in args:
                    args[args.index("keep")] = "all"
                elif "keep-and-load" in args:
                    args[args.index("keep-and-load")] = "upload"

            # fix for local directories: execute workflow in a shared
            # temporary directory. Note that restarting such a workflow
            # is not yet implemented.
            if self.workingdir_is_local:
                # change a default pipeline.log to remain in working dir as
                # otherwise the shared directory can not be cleaned up
                if options.logfile == "pipeline.log":
                    options.logfile = os.path.join(self.workingdir, options.logfile)

                self.shared_workingdir = P.get_temp_dir(shared=True)
                self.logger.debug("working directory is local, switching to shared".format(
                    self.shared_workingdir))
                self.logger.debug("working in shared temp dir {}".format(
                    self.shared_workingdir))
                os.chdir(self.shared_workingdir)
                P.PARAMS["workingdir"] = self.shared_workingdir

        return params

"""Exercise the helper against temporary Git repositories, without network access."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


HELPER = Path(__file__).resolve().parents[1] / 'scripts' / 'commit_and_push.sh'


class CommitAndPushTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='commit-push-test-', dir='/tmp')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.remote = self.root / 'remote.git'
        self.repo = self.root / 'checkout'
        self.repo.mkdir()
        self.git(self.root, 'init', '--bare', '--initial-branch=main', str(self.remote))
        self.git(self.repo, 'init', '--initial-branch=main')
        self.configure_author(self.repo)
        (self.repo / '.gitignore').write_text('.env\ncache/\n')
        shutil.copyfile(HELPER.parents[1] / '.gitattributes', self.repo / '.gitattributes')
        (self.repo / 'existing.txt').write_text('original\n')
        (self.repo / 'obsolete.txt').write_text('remove this file\n')
        self.git(self.repo, 'add', '-A')
        self.git(self.repo, 'commit', '-m', 'Initial fixture')
        self.git(self.repo, 'remote', 'add', 'origin', str(self.remote))
        self.git(self.repo, 'push', 'origin', 'main')
        (self.repo / 'scripts').mkdir()
        self.script = self.repo / 'scripts' / HELPER.name
        shutil.copyfile(HELPER, self.script)

    def git(self, directory, *args):
        result = subprocess.run(['git', '-C', str(directory), *args],
                                text=True, capture_output=True)
        if result.returncode:
            self.fail(f'git {args} failed:\n{result.stderr}')
        return result.stdout.strip()

    def configure_author(self, directory):
        self.git(directory, 'config', 'user.name', 'Automation integration test')
        self.git(directory, 'config', 'user.email', 'test@example.invalid')

    def helper(self, expect_success=True):
        # Invoke outside the checkout to exercise repository-root discovery.
        result = subprocess.run(['bash', str(self.script), 'Commit all analysis work'],
                                cwd=self.root, text=True, capture_output=True)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_new_modified_deleted_and_staged_files_are_pushed(self):
        (self.repo / 'existing.txt').write_text('modified\n')
        (self.repo / 'obsolete.txt').unlink()
        (self.repo / 'ora_analysis').mkdir()
        new_files = {
            'ora_analysis/analyze.py': 'print(51)\n',
            'ora_analysis/map.tsv': 'gene\tfunction\n',
            'ora_analysis/result.csv': 'gene,direction\n',
            'report.md': '# Analysis\n',
            'already_staged.txt': 'previously staged work\n',
        }
        for name, content in new_files.items():
            (self.repo / name).write_text(content)
        input_dir = self.repo / 'ora_analysis' / 'input'
        input_dir.mkdir()
        raw_input = b'gene,direction\r\nGene,up\r\n'
        (input_dir / 'data.csv').write_bytes(raw_input)
        self.git(self.repo, 'add', 'already_staged.txt')
        (self.repo / '.env').write_text('DUMMY_TEST_VALUE=not-a-credential\n')
        (self.repo / 'cache').mkdir()
        (self.repo / 'cache' / 'temp.txt').write_text('temporary\n')

        self.helper()
        pushed = set(self.git(self.remote, 'ls-tree', '-r', '--name-only', 'main').splitlines())
        self.assertTrue(set(new_files) <= pushed)
        self.assertIn('scripts/commit_and_push.sh', pushed)
        self.assertNotIn('obsolete.txt', pushed)
        self.assertNotIn('.env', pushed)
        self.assertNotIn('cache/temp.txt', pushed)
        self.assertEqual(self.git(self.remote, 'show', 'main:existing.txt'), 'modified')
        self.assertEqual(subprocess.check_output([
            'git', '-C', str(self.remote), 'show', 'main:ora_analysis/input/data.csv'
        ]), raw_input)
        self.assertEqual(self.git(self.repo, 'rev-parse', 'HEAD'), self.git(self.remote, 'rev-parse', 'main'))
        self.assertEqual(self.git(self.repo, 'status', '--porcelain'), '')

    def test_existing_commits_are_pushed_and_repeat_is_noop(self):
        self.git(self.repo, 'add', 'scripts/commit_and_push.sh')
        self.git(self.repo, 'commit', '-m', 'Unpushed local change')
        local = self.git(self.repo, 'rev-parse', 'HEAD')
        self.assertNotEqual(local, self.git(self.remote, 'rev-parse', 'main'))
        first = self.helper()
        self.assertIn('No new file changes', first.stdout)
        self.assertEqual(local, self.git(self.remote, 'rev-parse', 'main'))
        self.helper()
        self.assertEqual(local, self.git(self.repo, 'rev-parse', 'HEAD'))
        self.assertEqual(local, self.git(self.remote, 'rev-parse', 'main'))

    def test_diverged_remote_stops_before_staging_and_preserves_work(self):
        self.git(self.repo, 'add', 'scripts/commit_and_push.sh')
        self.git(self.repo, 'commit', '-m', 'Local-only commit')
        local = self.git(self.repo, 'rev-parse', 'HEAD')
        pending = self.repo / 'pending-analysis.md'
        pending.write_text('unsaved-to-Git analysis\n')
        (self.repo / 'existing.txt').write_text('working copy edits\n')
        self.git(self.repo, 'add', 'existing.txt')
        staged_before = self.git(self.repo, 'diff', '--cached')

        other = self.root / 'other-checkout'
        self.git(self.root, 'clone', str(self.remote), str(other))
        self.configure_author(other)
        (other / 'remote-update.md').write_text('another writer\n')
        self.git(other, 'add', '-A')
        self.git(other, 'commit', '-m', 'Remote advanced')
        self.git(other, 'push', 'origin', 'main')
        remote_before = self.git(self.remote, 'rev-parse', 'main')

        result = self.helper(expect_success=False)
        self.assertIn('commits missing locally', result.stderr)
        self.assertEqual(local, self.git(self.repo, 'rev-parse', 'HEAD'))
        self.assertEqual(remote_before, self.git(self.remote, 'rev-parse', 'main'))
        self.assertEqual(staged_before, self.git(self.repo, 'diff', '--cached'))
        self.assertEqual(pending.read_text(), 'unsaved-to-Git analysis\n')
        self.assertEqual((self.repo / 'existing.txt').read_text(), 'working copy edits\n')

    def test_real_whitespace_error_prevents_commit_and_push(self):
        local = self.git(self.repo, 'rev-parse', 'HEAD')
        remote = self.git(self.remote, 'rev-parse', 'main')
        bad = self.repo / 'bad.md'
        bad.write_text('unintended trailing space \n')
        result = self.helper(expect_success=False)
        self.assertIn('trailing whitespace', result.stdout)
        self.assertEqual(local, self.git(self.repo, 'rev-parse', 'HEAD'))
        self.assertEqual(remote, self.git(self.remote, 'rev-parse', 'main'))
        self.assertEqual(bad.read_text(), 'unintended trailing space \n')


if __name__ == '__main__':
    unittest.main()

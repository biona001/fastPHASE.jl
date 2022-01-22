using Tar, Inflate, SHA

download_urls = [
    "http://scheet.org/code/Linuxfp.tar.gz",
    "http://scheet.org/code/macOSfp.tar.gz"
    ]

for download_url in download_urls
    filename = download(download_url)
    println("sha256: ", bytes2hex(open(sha256, filename)))
    println("git-tree-sha1: ", Tar.tree_hash(IOBuffer(inflate_gzip(filename))))
end

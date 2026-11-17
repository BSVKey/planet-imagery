import struct, sys
def read_dbf(path):
    with open(path, 'rb') as f:
        hdr = f.read(32); n, hlen, rlen = struct.unpack('<xxxxIHH20x', hdr)
        fields = []
        while True:
            d = f.read(32)
            if d[0] == 0x0D: break
            fields.append((d[:11].split(b'\0')[0].decode(), chr(d[11]), d[16]))
        f.seek(hlen)
        for _ in range(n):
            rec = f.read(rlen)
            if rec[:1] == b'*': continue
            pos, row = 1, {}
            for name, t, ln in fields:
                row[name] = rec[pos:pos+ln].decode('utf-8', 'replace').strip(); pos += ln
            yield row
if __name__ == '__main__':
    rows = read_dbf(sys.argv[1]); r = next(rows); print(r)

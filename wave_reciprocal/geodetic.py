import numpy as np


def rotation_matrix(src_lat, src_lon):
    """ rotation matrix """
    theta = np.pi / 2. - np.radians(src_lat)
    phi = np.radians(src_lon)
    return np.array([
        [np.cos(theta) * np.cos(phi), -np.sin(phi), np.sin(theta) * np.cos(phi)],
        [np.cos(theta) * np.sin(phi), np.cos(phi), np.sin(theta) * np.sin(phi)],
        [-np.sin(theta), 0., np.cos(theta)]])


def lld_2_rtp(lld, r0):
    """ lat, lon, dep -> r, theta, phi """
    return np.array([r0 - lld[:, 2],
                     np.pi / 2. - np.radians(lld[:, 0]),
                     np.radians(lld[:, 1])]).T


def rtp_2_lld(rtp, r0):
    """ r, theta, phi -> lat, lon, dep """
    return np.array([np.degrees(np.pi / 2. - rtp[:, 1]),
                     np.degrees(rtp[:, 2]),
                     r0 - rtp[:, 0]]).T


def rtp_2_xyz(rtp):
    """ r, theta, phi -> x, y, z """
    x = rtp[:, 0] * np.sin(rtp[:, 1]) * np.cos(rtp[:, 2])
    y = rtp[:, 0] * np.sin(rtp[:, 1]) * np.sin(rtp[:, 2])
    z = rtp[:, 0] * np.cos(rtp[:, 1])
    return np.array([x, y, z]).T


def xyz_2_rtp(xyz, theta_where_undefined, phi_where_undefined):
    """ x, y, z -> r, theta, phi """
    # r
    r = np.linalg.norm(xyz, axis=1)
    # theta
    theta = np.arccos(xyz[:, 2] / r)
    theta[r == 0.] = theta_where_undefined[r == 0.]
    # phi
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    phi[phi < 0.] += np.pi
    phi[xyz[:, 1] < 0.] += np.pi
    phi[phi >= 2. * np.pi] -= 2. * np.pi
    s = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
    phi[s < 1e-8] = phi_where_undefined[s < 1e-8]
    return np.array([r, theta, phi]).T


class GeoPoints:
    """ GeoPoints class """
    r0 = 6371.

    def __init__(self, lld=np.zeros((1, 3)), phi_at_creation=None):
        self.lld = lld.copy()
        self.phi_at_creation = None
        if phi_at_creation is not None:
            self.phi_at_creation = phi_at_creation.copy()
        self.standardize_lon()

    def copy(self):
        return GeoPoints(self.lld, self.phi_at_creation)

    @staticmethod
    def create_rtp_src_centered(rtp_src, src_lat, src_lon):
        points = GeoPoints()
        points.set_rtp_src_centered(rtp_src, src_lat, src_lon)
        return points

    def standardize_lon(self, range0_360=True):
        lon = self.lld[:, 1]
        if range0_360:
            lon[lon < 0.] += 360.
        else:
            lon[lon > 180.] -= 360.

    ##############
    # geographic #
    ##############
    def get_rtp_geographic(self):
        return lld_2_rtp(self.lld, GeoPoints.r0)

    def set_rtp_geographic(self, rtp):
        self.lld = rtp_2_lld(rtp, GeoPoints.r0)
        self.standardize_lon()

    def get_xyz_geographic(self):
        return rtp_2_xyz(self.get_rtp_geographic())

    ###################
    # source-centered #
    ###################
    def get_xyz_src_centered(self, src_lat, src_lon):
        r_mat = rotation_matrix(src_lat, src_lon)
        xyz_geo = self.get_xyz_geographic()
        xyz_src = xyz_geo.dot(r_mat)
        return xyz_src

    def get_rtp_src_centered(self, src_lat, src_lon):
        xyz_src = self.get_xyz_src_centered(src_lat, src_lon)
        theta_where_undefined = np.zeros(len(self.lld))
        phi_where_undefined = self.phi_at_creation
        if phi_where_undefined is None:
            phi_where_undefined = np.zeros(len(self.lld))
        rtp_src = xyz_2_rtp(xyz_src, theta_where_undefined, phi_where_undefined)
        return rtp_src

    def set_rtp_src_centered(self, rtp_src, src_lat, src_lon):
        r_mat = rotation_matrix(src_lat, src_lon)
        xyz_src = rtp_2_xyz(rtp_src)
        xyz_geo = xyz_src.dot(r_mat.T)
        theta_where_undefined = np.zeros(len(rtp_src))
        self.phi_at_creation = rtp_src[:, 2].copy()
        rtp_geo = xyz_2_rtp(xyz_geo, theta_where_undefined, self.phi_at_creation)
        self.set_rtp_geographic(rtp_geo)

    def form_spz_frame(self, src_lat, src_lon, returns_phi_axis=False):
        rtp_me = self.get_rtp_src_centered(src_lat, src_lon)
        spz_me = np.array([rtp_me[:, 0] * np.sin(rtp_me[:, 1]),
                           rtp_me[:, 2],
                           rtp_me[:, 0] * np.cos(rtp_me[:, 1])]).T

        def spz_point(spz):
            r = np.sqrt(spz[:, 0] ** 2 + spz[:, 2] ** 2)
            rtp = np.array([r, np.arccos(spz[:, 2] / r), spz[:, 1]]).T
            return GeoPoints.create_rtp_src_centered(rtp, src_lat, src_lon)

        # ds, dz
        spz_ds = spz_me.copy()
        spz_dz = spz_me.copy()
        spz_ds[:, 0] += 1.
        spz_dz[:, 2] += 1.
        spz_ds_point = spz_point(spz_ds)
        spz_dz_point = spz_point(spz_dz)

        # frame
        vs = vector_xyz(self, spz_ds_point)
        vz = vector_xyz(self, spz_dz_point)
        vp = np.cross(vz, vs)
        if returns_phi_axis:
            return np.moveaxis(np.array([vs, vp, vz]), 0, 1), vp
        else:
            return np.moveaxis(np.array([vs, vp, vz]), 0, 1)

    def form_RTZ_frame(self, src_lat, src_lon):
        # vT
        _, vT = self.form_spz_frame(src_lat, src_lon, returns_phi_axis=True)

        # vZ
        rtp_me = self.get_rtp_src_centered(src_lat, src_lon)
        rtp_dr = rtp_me.copy()
        rtp_dr[:, 0] += 1.
        rtp_dr_point = GeoPoints.create_rtp_src_centered(rtp_dr, src_lat, src_lon)
        vZ = vector_xyz(self, rtp_dr_point)

        # vR
        vR = np.cross(vT, vZ)
        return np.moveaxis(np.array([vR, vT, vZ]), 0, 1)


#####################
# vector operations #
#####################
def vector_xyz(sp, ep, normalise=True):
    """ vector from sp to ep """
    v = ep.get_xyz_geographic() - sp.get_xyz_geographic()
    if normalise:
        v /= np.linalg.norm(v, axis=1)[:, None]
    return v


def vector_angle(v1, v2):
    """ angle between two vectors """
    inner = np.einsum('nd,nd->n', v1, v2)
    return np.arccos(inner / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)))


def rotate_vector_frame(vec_fr, x_fr, x_to):
    """ rotate vector """
    return vec_fr.dot(x_fr.dot(x_to.T))

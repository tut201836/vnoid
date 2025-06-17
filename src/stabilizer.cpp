#include "stabilizer.h"

#include "robot_base.h"
#include "rollpitchyaw.h"
#include "footstep.h"

namespace cnoid{
namespace vnoid{

Stabilizer::Stabilizer(){
    min_contact_force        = 1.0;
	force_ctrl_damping       = 0.0;
	force_ctrl_gain          = 0.0;
	force_ctrl_limit         = 0.0;
	moment_ctrl_damping      = 0.0;
	moment_ctrl_gain         = 0.0;
	moment_ctrl_limit        = 0.0;
	
	// default gain setting
	orientation_ctrl_gain_p = 10.0;
	orientation_ctrl_gain_d = 10.0;
	dcm_ctrl_gain           =  1.0;

	base_tilt_rate = 0.0;
	base_tilt_damping_p = 0.1;
	base_tilt_damping_d = 0.1;

	//
	recovery_moment_limit = 300.0;
	dcm_deviation_limit   = 0.3;
	
    for(int i = 0; i < 2; i++){
        dpos[i] = Vector3(0.0, 0.0, 0.0);
	    drot[i] = Vector3(0.0, 0.0, 0.0);
    }

}

void Stabilizer::CalcZmp(const Param& param, Centroid& centroid, vector<Foot>& foot) {
    const double eps = 1e-6;

    for (int i = 0; i < 2; i++) {
        const Vector3& f = foot[i].force;
        const Vector3& m = foot[i].moment;

        // 足裏のZ軸ベクトル（地面法線）
        Vector3 foot_z = foot[i].ori_ref.toRotationMatrix().col(2);

        // 法線方向の力（傾斜面に対応）
        double fz_local = foot_z.dot(f);
        double horiz_force_ratio = std::sqrt(f.x() * f.x() + f.y() * f.y()) / (fz_local + eps);

        // ヒステリシス付き接触判定 + 摩擦比による制約
        if (!foot[i].contact &&
            fz_local > contact_force_threshold_on &&
            horiz_force_ratio < friction_ratio_threshold) {
            foot[i].contact = true;
        } else if (foot[i].contact &&
                   (fz_local < contact_force_threshold_off ||
                    horiz_force_ratio > friction_ratio_threshold)) {
            foot[i].contact = false;
        }

        // スリップ検出（急激な法線力の減少）
        double prev_fz_local = foot[i].ori_ref.toRotationMatrix().col(2).dot(prev_force[i]);
        double fz_rate = fz_local - prev_fz_local;
        if (fz_rate < -30.0) {
            foot[i].contact = false;
        }
        prev_force[i] = f;

        // モーメントからZMP（ワールド座標）を計算
        Vector3 torque = m - foot[i].pos_ref.cross(f);
        Vector3 zmp_world = foot[i].pos_ref + foot_z.cross(torque) / (fz_local + eps);

        // ローカルZMPに変換
        if (foot[i].contact) {
            foot[i].zmp = foot[i].ori_ref.conjugate() * (zmp_world - foot[i].pos_ref);
        } else {
            foot[i].zmp = Vector3(0.0, 0.0, 0.0);
        }
    }

    // 両足とも非接地
    if (!foot[0].contact && !foot[1].contact) {
        foot[0].balance = foot[1].balance = 0.5;
        centroid.zmp = Vector3(0.0, 0.0, 0.0);
    } else {
        // 接地足の法線方向力に基づくZMP重み計算
        double fz0 = std::max(0.0, foot[0].ori_ref.toRotationMatrix().col(2).dot(foot[0].force));
        double fz1 = std::max(0.0, foot[1].ori_ref.toRotationMatrix().col(2).dot(foot[1].force));
        double fsum = fz0 + fz1 + eps;

        foot[0].balance = fz0 / fsum;
        foot[1].balance = fz1 / fsum;

        // ZMPをワールド座標で加重平均
        centroid.zmp =
            foot[0].balance * (foot[0].pos_ref + foot[0].ori_ref * foot[0].zmp) +
            foot[1].balance * (foot[1].pos_ref + foot[1].ori_ref * foot[1].zmp);
    }
}


void Stabilizer::CalcForceDistribution(const Param& param, Centroid& centroid, vector<Foot>& foot){
	// switch based on contact state
	if(!foot[0].contact_ref && !foot[1].contact_ref){
		foot[0].balance_ref = 0.5;
		foot[1].balance_ref = 0.5;
		foot[0].zmp_ref = Vector3(0.0, 0.0, 0.0);
		foot[1].zmp_ref = Vector3(0.0, 0.0, 0.0);
	}
	if( foot[0].contact_ref && !foot[1].contact_ref){
		foot[0].balance_ref = 1.0;
		foot[1].balance_ref = 0.0;
		foot[0].zmp_ref = foot[0].ori_ref.conjugate() * (centroid.zmp_ref - foot[0].pos_ref);
		foot[1].zmp_ref = Vector3(0.0, 0.0, 0.0);
	}
	if(!foot[0].contact_ref &&  foot[1].contact_ref){
		foot[0].balance_ref = 0.0;
		foot[1].balance_ref = 1.0;
		foot[0].zmp_ref = Vector3(0.0, 0.0, 0.0);
		foot[1].zmp_ref = foot[1].ori_ref.conjugate() * (centroid.zmp_ref - foot[1].pos_ref);
	}
	if( foot[0].contact_ref &&  foot[1].contact_ref){
		//
		Vector2 b;
		Vector3 pdiff  = foot[1].pos_ref - foot[0].pos_ref;
		double  pdiff2 = pdiff.squaredNorm();
		const double eps = 1.0e-10;
		if(pdiff2 < eps){
			b[0] = b[1] = 0.5;
		}
		else{
			b[0] = (pdiff.dot(foot[1].pos_ref - centroid.zmp_ref))/pdiff2;
			b[0] = std::min(std::max(0.0, b[0]), 1.0);
			b[1] = 1.0 - b[0];
		}

		foot[0].balance_ref = b[0];
		foot[1].balance_ref = b[1];

		Vector3 zmp_proj = b[0]*foot[0].pos_ref + b[1]*foot[1].pos_ref;

		double b2 = b.squaredNorm();
		foot[0].zmp_ref = (b[0]/b2) * (foot[0].ori_ref.conjugate() * (centroid.zmp_ref - zmp_proj));
		foot[1].zmp_ref = (b[1]/b2) * (foot[1].ori_ref.conjugate() * (centroid.zmp_ref - zmp_proj));
	}

	// limit zmp
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 3; j++){
			foot[i].zmp_ref[j] = std::min(std::max(param.zmp_min[j], foot[i].zmp_ref[j]), param.zmp_max[j]);
		}
	}

	for(int i = 0; i < 2; i++){
		// force and moment to realize desired Zmp
		foot[i].force_ref     =  foot[i].ori_ref.conjugate() * (foot[i].balance_ref * centroid.force_ref);
		foot[i].moment_ref[0] =  foot[i].force_ref.z() * foot[i].zmp_ref.y();
		foot[i].moment_ref[1] = -foot[i].force_ref.z() * foot[i].zmp_ref.x();
		foot[i].moment_ref[2] =  foot[i].balance_ref * centroid.moment_ref.z();
	}
}

void Stabilizer::CalcBaseTilt(const Timer& timer, const Param& param, Base& base, Vector3 theta, Vector3 omega){
	// desired angular acceleration for regulating orientation (in local coordinate)
	Vector3 omegadd_local(
		-(orientation_ctrl_gain_p*theta.x() + orientation_ctrl_gain_d*omega.x()),
		-(orientation_ctrl_gain_p*theta.y() + orientation_ctrl_gain_d*omega.y()),
		0.0
	);

	// 
	Vector3 omegadd_base(
		-base_tilt_rate*omegadd_local.x() - base_tilt_damping_p*base.angle_ref.x() - base_tilt_damping_d*base.angvel_ref.x(),
		-base_tilt_rate*omegadd_local.y() - base_tilt_damping_p*base.angle_ref.y() - base_tilt_damping_d*base.angvel_ref.y(),
		 0.0
	);

	base.angle_ref  += base.angvel_ref*timer.dt;
	base.ori_ref = FromRollPitchYaw(base.angle_ref);
	base.angvel_ref += omegadd_base*timer.dt;
}

void Stabilizer::CalcDcmDynamics(const Timer& timer, const Param& param, const Base& base, const vector<Foot>& foot, Vector3 theta, Vector3 omega, Centroid& centroid){
	double T = param.T;
	double m = param.total_mass;
	double h = param.com_height;
	double T_mh  = T/(m*h);
	double T2_mh = T*T_mh;

	Vector3 offset(0.0, 0.0, param.com_height);

	// desired angular acceleration for regulating orientation (in local coordinate)
	Vector3 omegadd_local(
		-(orientation_ctrl_gain_p*theta.x() + orientation_ctrl_gain_d*omega.x()),
		-(orientation_ctrl_gain_p*theta.y() + orientation_ctrl_gain_d*omega.y()),
		//-(orientation_ctrl_gain_p*theta.z() + orientation_ctrl_gain_d*omega.z())
		0.0
	);
	// desired moment (in local coordinate)
	Vector3 Ld_local(
		param.nominal_inertia.x()*omegadd_local.x(),
		param.nominal_inertia.y()*omegadd_local.y(),
		param.nominal_inertia.z()*omegadd_local.z()
	);
	// limit recovery moment for safety
	for(int i = 0; i < 3; i++){
		Ld_local[i] = std::min(std::max(-recovery_moment_limit, Ld_local[i]), recovery_moment_limit);
	}
	// desired moment (in global coordinate);
	Vector3 Ld = base.ori_ref * Ld_local;

	// virtual disturbance applied to DCM dynamics to generate desired recovery moment
	Vector3 delta = Vector3(-T_mh*Ld.y(), T_mh*Ld.x(), 0.0);

	// calc zmp for regulating dcm
	const double rate = 0.0;
	centroid.zmp_ref = centroid.zmp_target + dcm_ctrl_gain*(centroid.dcm_ref - centroid.dcm_target) + rate*T*delta;

	// project zmp inside support region
	if( ( foot[0].contact_ref && !foot[1].contact_ref) ||
		(!foot[0].contact_ref &&  foot[1].contact_ref) ){
		bool sup = (foot[0].contact_ref ? 0 : 1);
		Vector3 zmp_local = foot[sup].ori_ref.conjugate()*(centroid.zmp_ref - foot[sup].pos_ref);
		for(int j = 0; j < 3; j++){
			zmp_local[j] = std::min(std::max(param.zmp_min[j], zmp_local[j]), param.zmp_max[j]);
		}
		centroid.zmp_ref = foot[sup].pos_ref + foot[sup].ori_ref*zmp_local;
	}

	// calc DCM derivative
	Vector3 dcm_d = (1/T)*(centroid.dcm_ref - (centroid.zmp_ref + Vector3(0.0, 0.0, h))) + delta;

	// calc CoM acceleration
	centroid.com_acc_ref = (1/T)*(dcm_d - centroid.com_vel_ref);

	// update DCM
	centroid.dcm_ref += dcm_d*timer.dt;
	// limit deviation from reference dcm
	for(int j = 0; j < 3; j++){
		centroid.dcm_ref[j] = std::min(std::max(centroid.dcm_target[j] - dcm_deviation_limit, centroid.dcm_ref[j]), centroid.dcm_target[j] + dcm_deviation_limit);
	}

	// calc CoM velocity from dcm
	centroid.com_vel_ref = (1/T)*(centroid.dcm_ref - centroid.com_pos_ref);

	// update CoM position
	centroid.com_pos_ref += centroid.com_vel_ref*timer.dt;
}

void Stabilizer::InitializeState(Centroid& centroid, Base& base, std::vector<Foot>& foot, const Param& param) {
    double h = param.com_height;
    double T = param.T;

    // --- 足位置の平均でZMP/CoM初期化（両足接地仮定） ---
    Vector3 avg_foot_pos = 0.5 * (foot[0].pos_ref + foot[1].pos_ref);
    Vector3 com_pos = avg_foot_pos + Vector3(0.0, 0.0, h);

    centroid.com_pos_ref = com_pos;
    centroid.com_vel_ref = Vector3::Zero();
    centroid.dcm_ref     = com_pos;
    centroid.dcm_target  = centroid.dcm_ref;
    centroid.zmp_ref     = avg_foot_pos;
    centroid.zmp_target  = centroid.zmp_ref;

    // --- 足の状態初期化 ---
    for (int i = 0; i < 2; ++i) {
        foot[i].contact = true;
        foot[i].contact_ref = true;
        foot[i].force = Vector3(0.0, 0.0, param.total_mass * param.gravity / 2.0);
        foot[i].moment = Vector3::Zero();
        foot[i].force_ref = foot[i].force;
        foot[i].moment_ref = Vector3::Zero();
    }

    // --- Base姿勢初期化 ---
    base.angle = Vector3::Zero();
    base.angvel = Vector3::Zero();
    base.angle_ref = Vector3::Zero();
    base.angvel_ref = Vector3::Zero();
    base.ori_ref = FromRollPitchYaw(base.angle_ref);

    // --- ZMP補正変位初期化 ---
    for (int i = 0; i < 2; ++i) {
        dpos[i] = Vector3::Zero();
        drot[i] = Vector3::Zero();
        prev_force[i] = Vector3::Zero();
    }
}

void Stabilizer::Update(const Timer& timer, const Param& param, Centroid& centroid, Base& base, vector<Foot>& foot){
	// calc zmp from forces
    CalcZmp(param, centroid, foot);

	// error between desired and actual base link orientation
	Vector3 theta = base.angle  - base.angle_ref;
	Vector3 omega = base.angvel - base.angvel_ref;

	// calc base link tilt
	CalcBaseTilt(timer, param, base, theta, omega);
	
	// calc dcm and zmp 
	CalcDcmDynamics(timer, param, base, foot, theta, omega, centroid);

	// calc desired force applied to CoM
	// --- [重力ベクトル補正] 坂道対応 ---
	// ロボットのZ軸方向に合わせて重力を補正（滑りを防ぐ）
	Vector3 gravity_dir = base.ori_ref * Vector3(0.0, 0.0, -1.0);

	// NOTE: slope_gain を調整することで徐々に適応可能（0.0〜1.0）
	double slope_gain = 0.3;

	Vector3 gravity_vec_inclined = param.gravity * gravity_dir;
	Vector3 gravity_vec_flat     = Vector3(0.0, 0.0, param.gravity);
	Vector3 gravity_vector = slope_gain * gravity_vec_inclined + (1.0 - slope_gain) * gravity_vec_flat;

	centroid.force_ref = param.total_mass * (centroid.com_acc_ref + gravity_vector);

	centroid.moment_ref = Vector3(0.0, 0.0, 0.0);

	// calculate desired forces from desired zmp
	CalcForceDistribution(param, centroid, foot);

	for(int i = 0; i < 2; i++){
		// ground reaction force control
		if( foot[i].contact ){
			for(int j = 0; j < 3; j++){
				dpos[i][j] += (-force_ctrl_damping*dpos[i][j] + force_ctrl_gain*(foot[i].force_ref[j] - foot[i].force[j]))*timer.dt;
				dpos[i][j] = std::min(std::max(-force_ctrl_limit, dpos[i][j]), force_ctrl_limit);

				drot[i][j] += (-moment_ctrl_damping*drot[i][j] + moment_ctrl_gain*(foot[i].moment_ref[j] - foot[i].moment[j]))*timer.dt;
				drot[i][j] = std::min(std::max(-moment_ctrl_limit, drot[i][j]), moment_ctrl_limit);
			}

			// feedback to desired foot pose
			foot[i].pos_ref   += -dpos[i];
			foot[i].angle_ref += -drot[i];
      foot[i].ori_ref = FromRollPitchYaw(foot[i].angle_ref);
		}
	}

}


}
}
import hashlib
import hmac
import time

from pydantic import BaseModel

from .utils import api_client


class Organization(BaseModel):
    """
    A Rowan organization

    :ivar name: Name of the organization.
    :ivar weekly_credits: Weekly credits of the organization.
    :ivar credits: Credits of the organization.
    """

    name: str
    weekly_credits: float | None = None
    credits: float | None = None


class OrganizationRole(BaseModel):
    """
    A Rowan organization role

    :ivar name: Name of the organization role.
    """

    name: str


class SubscriptionPlan(BaseModel):
    """
    A Rowan subscription plan

    :ivar name: Name of the subscription plan.
    """

    name: str


class IndividualSubscription(BaseModel):
    """
    A Rowan individual subscription

    :ivar subscription_plan: Subscription plan of the individual subscription.
    """

    subscription_plan: SubscriptionPlan


class User(BaseModel):
    """
    A Rowan user

    :ivar uuid: UUID of the user.
    :ivar username: Username of the user.
    :ivar email: Email of the user.
    :ivar firstname: First name of the user.
    :ivar lastname: Last name of the user.
    :ivar weekly_credits: Weekly credits of the user.
    :ivar credits: Credits of the user.
    :ivar billing_name: Billing name of the user.
    :ivar billing_address: Billing address of the user.
    :ivar credit_balance_warning: Credit balance warning of the user.
    :ivar organization: Organization of the user.
    :ivar organization_role: Organization role of the user.
    :ivar individual_subscription: Individual subscription of the user.
    :ivar enabled_workflows: Workflow types this account may submit, as backend slugs that
        mostly match the ``submit_<slug>_workflow`` names (note ``molecular_dynamics`` is
        protein MD).
    :ivar feature_list: Feature flags enabled for this account.
    """

    uuid: str
    username: str
    email: str
    firstname: str | None = None
    lastname: str | None = None
    weekly_credits: float | None = None
    credits: float | None = None
    billing_name: str | None = None
    billing_address: str | dict | None = None
    credit_balance_warning: float | None = None
    organization: Organization | None = None
    organization_role: OrganizationRole | None = None
    individual_subscription: IndividualSubscription | None = None
    enabled_workflows: list[str] = []
    feature_list: list[str] = []

    @property
    def name(self) -> str:
        return f"{self.firstname or ''} {self.lastname or ''}".strip()

    @property
    def organization_name(self) -> str | None:
        return self.organization.name if self.organization else None

    @property
    def organization_weekly_credits(self) -> float | None:
        return self.organization.weekly_credits if self.organization else None

    @property
    def organization_credits(self) -> float | None:
        return self.organization.credits if self.organization else None

    @property
    def organization_role_name(self) -> str | None:
        return self.organization_role.name if self.organization_role else None

    @property
    def subscription_plan_name(self) -> str | None:
        return (
            self.individual_subscription.subscription_plan.name
            if self.individual_subscription
            else None
        )

    def credits_available_string(self) -> str:
        """
        Returns a string showing available credits, including organization credits if applicable

        :returns: String showing available credits
        """
        individual_credits = f"Weekly Credits: {self.weekly_credits}\nCredits: {self.credits}"
        if self.organization is not None:
            return (
                individual_credits
                + f"\nOrganization Weekly Credits: {self.organization_weekly_credits}\n"
                f"Organization Credits: {self.organization_credits}"
            )

        return individual_credits


def whoami() -> User:
    """
    Returns the current user
    """

    with api_client() as client:
        response = client.get("/user/me")
        response.raise_for_status()

    return User(**response.json())


def get_webhook_secret() -> str | None:
    """Get current webhook secret, or None if not set."""
    with api_client() as client:
        response = client.get("/user/me/webhook_secret")
        response.raise_for_status()
    return response.json().get("webhook_secret")


def create_webhook_secret() -> str:
    """Create a webhook secret if one doesn't exist, return it."""
    with api_client() as client:
        response = client.post("/user/me/webhook_secret")
        response.raise_for_status()
    return response.json()["webhook_secret"]


def rotate_webhook_secret() -> str:
    """Generate a new webhook secret, invalidating the old one."""
    with api_client() as client:
        response = client.post("/user/me/webhook_secret/rotate")
        response.raise_for_status()
    return response.json()["webhook_secret"]


def verify_webhook_secret(
    raw_body: bytes,
    signature_header: str,
    secret: str,
    max_age_seconds: int = 300,
) -> bool:
    """Verify an incoming webhook request from Rowan.

    :param raw_body: Raw (unparsed) request body bytes.
    :param signature_header: Value of the X-Rowan-Signature header.
    :param secret: Your webhook secret (from :func:`create_webhook_secret` or
        :func:`rotate_webhook_secret`).
    :param max_age_seconds: Reject requests older than this many seconds (default 5 min).
    :returns: True if the signature is valid and the request is fresh.
    """
    try:
        parts = dict(part.split("=", 1) for part in signature_header.split(","))
    except ValueError:
        return False
    timestamp = parts.get("t", "")
    received_digest = parts.get("sha256", "")

    if not timestamp or not received_digest:
        return False

    if abs(time.time() - int(timestamp)) > max_age_seconds:
        return False

    expected_digest = hmac.new(
        secret.encode(),
        f"{timestamp}.".encode() + raw_body,
        hashlib.sha256,
    ).hexdigest()

    return hmac.compare_digest(expected_digest, received_digest)
